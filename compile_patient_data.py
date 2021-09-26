from typing import IO, Iterable, Tuple
import pandas as pd
import os
import csv
import re
import json

json.dumps_nws = lambda obj: json.dumps(obj, separators=(',', ':'))

FILE_NAME_PATTERN = re.compile(r"(\w{4}-\w{2}-\w{4})-(.+)\.(.+)+\.txt")

DUPE_miRNA_PATTERN = re.compile(r"(?i)(hsa-\w+-\d+[a-z]?)-\d+")
BASE_miRNA_PATTERN = re.compile(r"(?i)(hsa-\w+-\d+[a-z]?)")

GENCODE_FILE = 'gencode.gene.info.v22.tsv'
GENE_TYPES_OF_INTEREST = {'protein_coding', }
                          # 'snRNA', 'rRNA',  # 'miRNA',
                          # 'snoRNA', 'ribozyme', 'macro_lncRNA', }  # 'vaultRNA', 'Mt_tRNA', 'Mt_rRNA'}

MIRBASE_MATURE_OLD_FILE = 'mature.v17.fa.gz'
MIRBASE_MATURE_NEW_FILE = 'mature.v21.fa.gz'


def source_sample(ancillary_id: str = None, is_normal: bool = False):
    if ancillary_id:
        if ancillary_id.startswith('01'):
            return 'PRIMARY_TUMOR'
        elif ancillary_id.startswith('11'):
            return 'SOLID_TISSUE_NORMAL'
    elif is_normal:
        return 'SOLID_TISSUE_NORMAL'
    else:
        return 'PRIMARY_TUMOR'


def get_seq_file_params(filename):
    m = FILE_NAME_PATTERN.match(filename)
    if m:
        patient_id = m.group(1)
        ancillary_id = m.group(2)
        # entity_id = f'{patient_id}-{ancillary_id}'
        seq_type, count_type = m.group(3).split('-')
        source = source_sample(ancillary_id)
        return (patient_id, seq_type, source)


def unpaired_seq_files(data_dir):
    sample_set = {}
    verified_sample_set = set()

    verified_sample_code = ('miRNA', 'mRNA')

    for member in os.listdir(data_dir):
        seq_file_params = get_seq_file_params(member)
        if seq_file_params:
            (patient_id, seq_type, source) = seq_file_params
            print(f"Processing file: {patient_id}, {seq_type}")

            if patient_id in sample_set:
                sample_set[patient_id] += (seq_type,)
            else:
                sample_set[patient_id] = (seq_type,)
            if sample_set[patient_id] == verified_sample_code:
                verified_sample_set.add(patient_id)
    total = len(sample_set)
    count = 0
    unpaired_samples = set()
    for sample in sample_set.keys():
        if sample not in verified_sample_set:
            unpaired_samples.add(sample)
            print(f"{sample} has only {sample_set[sample]} seq data")
        else:
            count += 1
    print(f"Verified {count} of {total} samples")
    return unpaired_samples


def init_gencode_name_map(gencode_map_file: str = GENCODE_FILE, strip_version: bool = False):
    # Initialize GenCode to Gene Name map
    with open(gencode_map_file, 'r', newline='') as map_file:
        map_csv = csv.reader(map_file, dialect='excel-tab')
        next(map_csv)  # skip header
        ENSG_CODE_TO_NAME_MAP = {}
        for row in map_csv:
            # Format (TSV):
            #    0          1          2      3      4       5        6           7            8           9           10        11
            # gene_id   gene_name   seqname start   end   strand  gene_type   gene_status havana_gene full_length exon_length exon_num
            if row[7] == 'KNOWN' and row[6] in GENE_TYPES_OF_INTEREST:
                gene_id = row[0]
                gene_id = gene_id[:gene_id.index('.')] if strip_version else gene_id
                ENSG_CODE_TO_NAME_MAP[gene_id] = row[1]
        return ENSG_CODE_TO_NAME_MAP


def read_mirbase_fasta(file):
    data = []
    line_buffer = []

    import pathlib
    import gzip

    file_ext = pathlib.Path(file).suffix.lower()
    is_gzip = file_ext in ('.gz', '.gzip')
    with (gzip.open(file, 'rt') if is_gzip else open(file, 'r')) as fa_file:
        for line in fa_file:
            if line.startswith(">"):
                line = line[1:].strip().split()
                # >name accession genus species name_w/o_species_code
                line = [line[0], line[1], f"{line[2]} {line[3]}"]
                line_buffer.extend(line)
            else:
                # sequence
                line_buffer.append(line.strip())
                data.append(line_buffer)
                line_buffer = []
    data_df = pd.DataFrame(data, columns=('Name', 'Accession', 'Species', 'Sequence'))
    data_df.set_index('Accession', inplace=True)
    return data_df


def init_mirbase_name_map(mirbase_mature_old_file: str = MIRBASE_MATURE_OLD_FILE,
                          mirbase_mature_new_file: str = MIRBASE_MATURE_NEW_FILE,
                          filter_species: str = None, preferred_name: str = 'Old Name'):
    mirbase_old = read_mirbase_fasta(mirbase_mature_old_file)

    if filter_species:
        mirbase_old = mirbase_old[mirbase_old['Species'] == filter_species]

    mirbase_new = read_mirbase_fasta(mirbase_mature_new_file)

    if filter_species:
        mirbase_new = mirbase_new[mirbase_new['Species'] == filter_species]

    mirbase_df = mirbase_new
    mirbase_df.rename({'Name': 'New Name'}, axis='columns', inplace=True)
    mirbase_df['Old Name'] = mirbase_df.index.map(lambda acc: mirbase_old['Name'].loc[acc] if acc in mirbase_old.index else mirbase_df['New Name'].loc[acc])

    mirbase_dict = mirbase_df[preferred_name].to_dict()
    return mirbase_dict


class PatientDataAggregator():
    def __init__(self, gene_name_map_file=GENCODE_FILE,
                 mirbase_mature_old_file: str = MIRBASE_MATURE_OLD_FILE,
                 mirbase_mature_new_file: str = MIRBASE_MATURE_NEW_FILE,
                 preferred_name: str = 'Old Name') -> None:
        self.gene_name_map = init_gencode_name_map(gene_name_map_file)
        self.mirna_name_map = init_mirbase_name_map(mirbase_mature_old_file=mirbase_mature_old_file,
                                                    mirbase_mature_new_file=mirbase_mature_new_file,
                                                    filter_species='Homo sapiens',
                                                    preferred_name=preferred_name)

    def patient_data_as_df_row(self, patient_id: str, source: str, is_mutant: bool,
                               mrna_fpkm_file: IO, mirna_counts_file: IO) -> Tuple[pd.DataFrame, pd.DataFrame]:
        if (not mrna_fpkm_file) or (not mirna_counts_file):
            return None, None

        mrna_fpkm_data: pd.DataFrame = pd.read_csv(
            mrna_fpkm_file, sep='\t', index_col=False, header=0, usecols=[0, 1], names=('Gene', 'FPKM'))

        mrna_fpkm_data['Gene'] = mrna_fpkm_data['Gene'].apply(
            lambda ensg_id: self.gene_name_map.get(ensg_id))
        mrna_fpkm_data = mrna_fpkm_data.dropna()
        mrna_fpkm_data.set_index('Gene', inplace=True)

        patient_mrna_data = mrna_fpkm_data.transpose()
        patient_mrna_data.insert(loc=0, column='Patient_ID', value=[patient_id])
        patient_mrna_data.insert(loc=1, column='Source', value=[source])
        patient_mrna_data.insert(loc=2, column='Is_Mutant', value=[is_mutant])
        patient_mrna_data.set_index('Patient_ID', inplace=True)

        mirna_rpm_data: pd.DataFrame = pd.read_csv(
            mirna_counts_file, sep='\t', index_col=False, header=0, usecols=[0, 3, 5], names=('miRNA', 'RPM', 'miRNA_region'))

        def get_accession(mirna_region):
            if ',' in mirna_region:
                return mirna_region.split(',')[1]
            else:
                return None

        mirna_rpm_data['accession'] = mirna_rpm_data['miRNA_region'].apply(get_accession)
        mirna_rpm_data = mirna_rpm_data[mirna_rpm_data['accession'].apply(lambda x: not pd.isna(x))]

        mirna_rpm_data = mirna_rpm_data.dropna()
        mirna_rpm_data_deduped = mirna_rpm_data.groupby('accession')

        patient_mirna_data = {
            'Patient_ID': [patient_id],
            'Source': [source],
            'Is_Mutant': [is_mutant],
        }
        for key in mirna_rpm_data_deduped.groups:
            group = mirna_rpm_data_deduped.get_group(key)
            patient_mirna_data[self.mirna_name_map.get(group['accession'].iloc[0], group['miRNA'].iloc[0])] = group['RPM'].sum()

        patient_mirna_data = pd.DataFrame(patient_mirna_data)
        patient_mirna_data.set_index('Patient_ID', inplace=True)

        return patient_mrna_data, patient_mirna_data

    @staticmethod
    def merge_patient_data_rows(patient_data_rows: Iterable[Tuple[pd.DataFrame, pd.DataFrame]]) -> Tuple[pd.DataFrame, pd.DataFrame]:
        mrna_data_rows, mirna_data_rows = zip(*patient_data_rows)
        aggregate_mrna_data = pd.concat(mrna_data_rows).fillna(0.0)
        aggregate_mirna_data = pd.concat(mirna_data_rows).fillna(0.0)
        return aggregate_mrna_data, aggregate_mirna_data


def main(data_dir, output_file_name=None):

    ENSG_CODE_TO_NAME_MAP = init_gencode_name_map()

    gene_miRNA_data = {}
    entity_names = set()
    metadata = set()

    if not data_dir:
        data_dir = './'

    for member in os.listdir(data_dir):
        seq_file_params = get_seq_file_params(member)
        if seq_file_params:
            (patient_id, seq_type, source) = seq_file_params
            print(f"Processing file: {patient_id}, {seq_type}")
            metadata.add((patient_id, source))

            member_path = os.path.join(data_dir, member)
            with open(member_path, 'r', newline='') as input_file:
                input_csv = csv.reader(input_file, dialect='excel-tab')
                next(input_csv)  # skip header
                entity_name_value_pairs = {}
                if seq_type == 'miRNA':
                    for row in input_csv:
                        # rpm = reads per million
                        name, rpm = row[0], float(row[2])
                        entity_name_value_pairs[name] = rpm
                elif seq_type == 'mRNA':
                    for row in input_csv:
                        name, fpkm = row[0], float(row[1])
                        # get gene name corresponding to ENSG ID if relevant or skip
                        if name in ENSG_CODE_TO_NAME_MAP:
                            name = ENSG_CODE_TO_NAME_MAP[name]
                            entity_name_value_pairs[name] = fpkm

                gene_miRNA_data_key = (patient_id, source)
                entity_names.update(entity_name_value_pairs.keys())

                if gene_miRNA_data_key not in gene_miRNA_data:
                    gene_miRNA_data[gene_miRNA_data_key] = entity_name_value_pairs
                else:
                    gene_miRNA_data[gene_miRNA_data_key].update(entity_name_value_pairs)

    for key, entity_names_values in gene_miRNA_data.items():
        current_entity_names = entity_names_values.keys()
        if len(current_entity_names) < len(entity_names):
            for entity_name in entity_names - current_entity_names:
                # pad out missing values
                entity_names_values[entity_name] = 0.0
        entity_names_values = list(entity_names_values.items())
        entity_names_values.sort(key=lambda x: x[0])
        gene_miRNA_data[key] = entity_names_values

    print(f"Processed {len(metadata)} unique files")
    output_file_name = output_file_name or "compiled_data.csv"
    output_file_path = os.path.join(data_dir, output_file_name)
    with open(output_file_path, 'w', newline='') as output_file:
        output_csv = csv.writer(output_file)
        # Header
        name_list = list(entity_names)
        name_list.sort()
        meta_info = ['Patient ID', 'Source']
        header = meta_info + name_list
        count, total = 0, len(gene_miRNA_data)
        output_csv.writerow(header)
        for key in gene_miRNA_data:
            row = [key[0], key[1]]
            for entity_name, value in gene_miRNA_data[key]:
                row.append(value)
            count += 1
            print(f"Wrote {len(row)} entities for row {count} of {total}")
            output_csv.writerow(row)


if __name__ == '__main__':
    import sys
    main(sys.argv[1] if len(sys.argv) > 1 else None)
