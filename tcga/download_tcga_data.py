from typing import Tuple, List
import requests
import urllib3
import json
import gzip
import os
import io
import multiprocessing as mp

from gdc_explore_goi_mut import goi_mutants_file
from compile_patient_data import PatientDataAggregator, source_sample

DEFAULT_OUTPUT_DIR = 'data'

cases_endpt = "https://api.gdc.cancer.gov/cases"
files_endpt = "https://api.gdc.cancer.gov/files"
data_endpt = "https://api.gdc.cancer.gov/data"

normal_tissue_data = False

datasets = {
    'luad': {
        'program_name': "TCGA",
        'project_id': "TCGA-LUAD",
        'disease_type': "adenomas and adenocarcinomas",
        'primary_site': "bronchus and lung",
    },
    'gbm': {
        'program_name': "TCGA",
        'disease_type': "gliomas",
        'primary_site': "brain",
    },
    'hpcc': {
        'program_name': "TCGA",
        'disease_type': "adenomas and adenocarcinomas",
        'primary_site': "liver and intrahepatic bile ducts",
    },
    'brca': {
        'program_name': "TCGA",
        'project_id': "TCGA-BRCA",
        'disease_type': "ductal and lobular neoplasms",
        'primary_site': "breast",
    },
}


def get_case_filters(dataset_specifiers: dict, is_normal: bool = False) -> dict:
    base_params = {"op": "and",
                   "content": [
                       {
                           "op": "in",
                           "content": {
                               "field": "cases.samples.sample_type",
                               "value": [
                                   "solid tissue normal" if is_normal else "primary tumor"
                               ]
                           }
                       },
                       {
                           "op": "in",
                           "content": {
                               "field": "cases.summary.experimental_strategies.experimental_strategy",
                               "value": [
                                   "RNA-Seq",
                                   "miRNA-Seq"
                               ]
                           }
                       }
                   ]
                   }
    if program_name := dataset_specifiers.get('program_name'):
        base_params['content'].append({
                                      "op": "in",
                                      "content": {
                                          "field": "cases.project.program.name",
                                          "value": [program_name]
                                      }})
    if project_id := dataset_specifiers.get('project_id'):
        base_params['content'].append({
                                      "op": "in",
                                      "content": {
                                          "field": "cases.project.project_id",
                                          "value": [project_id]
                                      }})
    if disease_type := dataset_specifiers.get('disease_type'):
        base_params['content'].append({
                                      "op": "in",
                                      "content": {
                                          "field": "cases.disease_type",
                                          "value": [disease_type]
                                      }})
    if primary_site := dataset_specifiers.get('primary_site'):
        base_params['content'].append({
                                      "op": "in",
                                      "content": {
                                          "field": "cases.primary_site",
                                          "value": [primary_site]
                                      }})
    return base_params


def get_file_filters(case_id, is_normal=False):
    return {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {
                    "field": "cases.case_id",
                    "value": [
                        str(case_id)
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "cases.samples.sample_type",
                    "value": [
                        "solid tissue normal" if is_normal else "primary tumor"
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.data_type",
                    "value": [
                        "Gene Expression Quantification",
                        "Isoform Expression Quantification"
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.analysis.workflow_type",
                    "value": [
                        "HTSeq - FPKM",
                        "BCGSC miRNA Profiling"
                    ]
                }
            }
        ]
    }


CASE_ID_CODE = 'submitter_id'
ENTITY_ID_CODE = 'submitter_id'
FILE_NAME_CODE = 'file_name'


def get_hits(response):
    data = json.loads(response.content.decode("utf-8"))
    if data.get('warnings'):
        print(data['warnings'])
    return data["data"]["hits"]


def make_params(filters, *fields, size_lim=5000):
    return {
        "filters": json.dumps(filters),
        "fields": ','.join(fields),
        "format": "JSON",
        "size": str(size_lim)
    }


def paired_seq_file_entries(file_entries):
    if len(file_entries) < 2:
        return None
    paired_files = {
        'miRNA': [
            file_entry for file_entry in file_entries
            if '.isoforms.' in file_entry[FILE_NAME_CODE].lower()
        ],
        'mRNA': [
            file_entry for file_entry in file_entries
            if '.fpkm.' in file_entry[FILE_NAME_CODE].lower()
        ]
    }

    if paired_files['miRNA'] and paired_files['mRNA']:
        paired_files['miRNA'] = paired_files['miRNA'][0]
        paired_files['mRNA'] = paired_files['mRNA'][0]
        return paired_files
    else:
        return None


def process_files_for_case_entry(params: Tuple[int, PatientDataAggregator, List[str], dict]):
    num, patient_data_aggregator, mutant_case_ids, case_entry = params
    case_uuid = case_entry['id']
    case_id = case_entry[CASE_ID_CODE]

    file_params = make_params(get_file_filters(
        case_uuid), ENTITY_ID_CODE, FILE_NAME_CODE)

    print(f"Requesting files for Case ID {case_id}, No. {num}")
    files_response = requests.get(files_endpt, params=file_params)

    file_entries = get_hits(files_response)

    paired_file_entries = paired_seq_file_entries(file_entries)

    if not paired_file_entries:
        return (None, None)

    def get_file_data(file_entry):
        file_id = file_entry['id']
        full_file_name = file_entry[FILE_NAME_CODE]
        _, ext = os.path.splitext(full_file_name)
        params = json.dumps({"ids": file_id})
        decompress = ext.lower().endswith('.gz')

        print(f"Downloading file {full_file_name} for Case ID {case_id}, No. {num}")
        try:
            response = requests.post(data_endpt, data=params, headers={
                                     "Content-Type": "application/json"})
            file_data = io.BytesIO(response.content)
            file_data = gzip.GzipFile(fileobj=file_data) if decompress else file_data
            # print(f"Downloaded file {full_file_name} for Case ID {case_id}, No. {num}")
            return file_data
        except (requests.exceptions.ConnectionError,
                urllib3.exceptions.NewConnectionError,
                urllib3.exceptions.MaxRetryError,
                OSError):
            print(f"Failed to download file {full_file_name} for Case ID {case_id}, No. {num}; skipping.")
            return None

    paired_file_entries = {key: get_file_data(value)
                           for key, value in paired_file_entries.items()}

    patient_data = patient_data_aggregator.patient_data_as_df_row(patient_id=case_id,
                                                                  source=source_sample(is_normal=normal_tissue_data),
                                                                  is_mutant=(case_id in mutant_case_ids),
                                                                  mrna_fpkm_file=paired_file_entries['mRNA'],
                                                                  mirna_counts_file=paired_file_entries['miRNA'])

    print(f"Processed data for Case ID {case_id}, No. {num}")
    return patient_data


def main(dataset_name, output_dir):
    if not output_dir:
        output_dir = f"{dataset_name}_{DEFAULT_OUTPUT_DIR}".upper()

    patient_data_aggregator = PatientDataAggregator()
    case_filters = get_case_filters(datasets[dataset_name], normal_tissue_data)
    case_params = make_params(case_filters, CASE_ID_CODE)
    print(f"Requesting entries from TCGA for dataset {dataset_name}")
    cases_response = requests.get(cases_endpt, params=case_params)
    case_hits = get_hits(cases_response)
    total_cases = len(case_hits)
    print(f"Total {dataset_name} cases = {total_cases}")

    # print(f"Loading list of entries with mutations in genes of interest from file {goi_mutants_file}")
    mutant_case_ids = None
    with open(goi_mutants_file, 'r') as mutant_cases:
        mutant_case_ids = [mutant[CASE_ID_CODE] for mutant in json.load(mutant_cases)]

    print(f"Loaded {len(mutant_case_ids)} entries with deleterious mutations in genes of interest from file {goi_mutants_file}")

    # non_mutant_case_hits = [case_hit for case_hit in case_hits if case_hit[CASE_ID_CODE] not in goi_mutant_case_ids]

    # case_hits = non_mutant_case_hits

    # print(f"Final {dataset_name} cases = {len(non_mutant_case_hits)}")

    with mp.Pool() as worker_pool:
        def on_complete(result):
            aggregate_mrna_data, aggregate_mirna_data = patient_data_aggregator.merge_patient_data_rows(result)
            aggregate_mrna_data.to_csv(f'{output_dir}/{dataset_name}_mRNA.csv.gz', compression='gzip')
            aggregate_mirna_data.to_csv(f'{output_dir}/{dataset_name}_miRNA.csv.gz', compression='gzip')
            print("Data Download Complete, Written Processed CSV Files for Matched RNA-Seq Data")

        force_abort = False

        def on_error(error):
            nonlocal force_abort
            force_abort = True
            raise error

        cases = len(case_hits)
        completion = worker_pool.map_async(process_files_for_case_entry,
                                           zip(range(1, len(case_hits) + 1),
                                               [patient_data_aggregator] * cases,
                                               [mutant_case_ids] * cases, case_hits),
                                           callback=on_complete,
                                           error_callback=on_error)
        while not completion.ready():
            if force_abort:
                sys.exit(-1)
            completion.wait(1)


if __name__ == '__main__':
    import sys
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
