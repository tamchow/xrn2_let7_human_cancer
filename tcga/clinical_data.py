#%%
from pathlib import Path
import numpy as np
import requests
import json
import pandas as pd

def process_project_clinical_data(patient_info_file):
    basis_col = "XRN2"
    patient_info = pd.read_csv(patient_info_file, usecols=["Patient_ID", basis_col])
    basis_median = patient_info[basis_col].median()
    basis_group_col = f"{basis_col} Expr Group"
    patient_info[basis_group_col] = patient_info[basis_col].apply(
        lambda expr: "low"
        if expr < basis_median
        else "high"
        if expr > basis_median
        else "normal"
    )

    GDC_CASES_ENDPNT = "https://api.gdc.cancer.gov/cases"


    def get_clinical_data(submitter_ids):
        submitter_ids = list(set(submitter_ids))
        filter = {"op": "in", "content": {"field": "submitter_id", "value": submitter_ids}}
        params = {
            "filters": json.dumps(filter),
            "expand": "diagnoses,demographic",
            "size": len(submitter_ids),
        }
        response = requests.get(GDC_CASES_ENDPNT, params=params)
        hits = response.json()
        clinical_data = {}
        try:
            hits = hits["data"]["hits"]
            for hit in hits:
                if "diagnoses" not in hit or "demographic" not in hit:
                    print(f'Patient ID {hit["submitter_id"]} has no diagnoses or demographic')
                diagnoses = hit.get("diagnoses", None)
                demographic = hit.get("demographic", None)
                clinical_data[hit["submitter_id"]] = diagnoses, demographic
            return clinical_data
        except KeyError as keyError:
            print(keyError.args)
            print(json.dumps(hits, indent=2))


    CLINICAL_DATA_COLS = [
        "age",
        "primary_diagnosis",
        # "ajcc_pathologic_stage",
        "prior_malignancy",
        "synchronous_malignancy",
        "prior_treatment",
        # "progression_or_recurrence",
        "race",
        "sex",
        "alive",
    ]


    def ifna(val: float) -> float:
        return np.nan if pd.isna(val) else val


    def interesting_clinical_data(diagnoses, demographic):
        age = np.nan
        primary_diagnosis = pd.NA
        # ajcc_pathologic_stage = pd.NA
        prior_malignancy = pd.NA
        synchronous_malignancy = pd.NA
        prior_treatment = pd.NA
        # progression_or_recurrence = pd.NA
        if diagnoses:
            diagnosis = diagnoses[0]
            age = ifna(diagnosis["age_at_diagnosis"]) / 365  # convert from days to years
            primary_diagnosis = diagnosis["primary_diagnosis"].lower()
            # ajcc_pathologic_stage = [
            #     diagnosis["ajcc_pathologic_stage"].lower() for diagnosis in diagnoses
            # ]
            prior_malignancy = diagnosis["prior_malignancy"].lower() == "yes"
            synchronous_malignancy = diagnosis["synchronous_malignancy"].lower() == "yes"
            prior_treatment = diagnosis["prior_treatment"].lower() == "yes"
            # progression_or_recurrence = [
            #     diagnosis["progression_or_recurrence"] for diagnosis in diagnoses
            # ]
        
        race = pd.NA
        sex = pd.NA
        alive = pd.NA
        if demographic:
            race = demographic.get("race", "not reported").lower()
            sex = demographic.get("gender", "not reported").lower()
            alive = demographic.get("vital_status", "not reported").lower() == "alive"
        clinical_data = [
            age,
            primary_diagnosis,
            # ajcc_pathologic_stage,
            prior_malignancy,
            synchronous_malignancy,
            prior_treatment,
            # progression_or_recurrence,
            race,
            sex,
            alive,
        ]
        return pd.Series(clinical_data, index=CLINICAL_DATA_COLS)


    all_clinical_data = get_clinical_data(patient_info["Patient_ID"])

    patient_info = patient_info.join(
        patient_info["Patient_ID"].apply(
            lambda submitter_id: interesting_clinical_data(*all_clinical_data[submitter_id])
        )
    )

    for col in patient_info.columns:
        if patient_info[col].dtype == "object":
            patient_info[col] = patient_info[col].astype("category")

    with pd.ExcelWriter(
        f"{patient_info_file.parent}/{patient_info_file.stem}.xlsx"
    ) as out_excel:
        patient_info.sort_values(by=[basis_group_col], inplace=True)
        patient_info.set_index('Patient_ID', inplace=True)
        patient_info.to_excel(
            out_excel, sheet_name="clinical_data", index_label='Submitter ID'
        )

        for group, df in patient_info.groupby(basis_group_col):
            df.describe(include="all").to_excel(
                out_excel, sheet_name=f"{basis_col} {group}"
            )

projects = ['LUAD', 'GBM', 'HPCC']
for project in projects:
    process_project_clinical_data(Path(f"{project}_DATA/{project.lower()}_sample_subset.csv.gz"))
#%%
