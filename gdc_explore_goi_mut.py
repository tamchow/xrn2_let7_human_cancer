import urllib.parse as up
import json

from compile_patient_data import init_gencode_name_map

goi_mutants_file = "goi_mutants.json"

all_genes_of_interest = ['XRN2', 'CDKN2AIP', 'NKRF', 'NOP56', 'NOP58',
                     'XRN1', 'DCPS', 'DCP1A', 'DCP1B', 'DCP2', 'SND1', 'ZSWIM8']

# vep_impacts = ["moderate", "high"]
sift_impacts = ["deleterious"]
polyphen_impacts = ["probably_damaging", "possibly_damaging"]


def main(genes_of_interest):
    gene_name_to_ensg_id_map = {name: ensg_id for ensg_id, name in init_gencode_name_map(strip_version=True).items()}

    explore_query = {
        "op": "and",
        "content": [
        {
            "op": "in",
            "content": {
                "field": "cases.samples.sample_type",
                "value": ["primary tumor"]
            }
        },
        {
            "op": "in",
            "content": {
                "field": "genes.gene_id",
                "value": [gene_name_to_ensg_id_map[gene] for gene in genes_of_interest]
            }
        },
        # {
        #     "op": "in",
        #     "content": {
        #         "field": "ssms.consequence.transcript.annotation.vep_impact",
        #         "value": vep_impacts
        #     }
        # },
        {
            "op": "in",
            "content": {
                "field": "ssms.consequence.transcript.annotation.sift_impact",
                "value": sift_impacts
            }
        },
        {
            "op": "in",
            "content": {
                "field": "ssms.consequence.transcript.annotation.polyphen_impact",
                "value": polyphen_impacts
            }
        }]
    }

    SAFE_QUOTE_CHARS = '{\"}'
    gdc_explore_url = f"https://portal.gdc.cancer.gov/exploration?filters={up.quote_plus(json.dumps_nws(explore_query), safe=SAFE_QUOTE_CHARS)}"

    print("Please copy and paste the following URL and press the 'JSON' button to download the list of patients having mutations in genes of interest.")
    print(f"Please save the downloaded JSON file in the same directory as all these scripts with the name '{goi_mutants_file}'")
    print(gdc_explore_url)


if __name__ == '__main__':
    import sys
    main(sys.argv[1:] if len(sys.argv) > 1 else all_genes_of_interest)
