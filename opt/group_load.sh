#!/usr/bin/env bash

if [[ $# -ne 1 ]]
    then
        echo "usage: $0 config_file"
        echo "  config_file: path to the file containing environment variables for the patient group to be used"
        exit 2
fi

if [[ ! -f $1 ]]
    then
        echo "invalid file path: $1"
        exit 2
fi

source "$1"

list_idx=0
for i in {00001..00300}
    do
        candig_repo add-variantset -r \
            -I "${GROUP}"/"${GROUP}"_${i}.vcf.gz.tbi \
            -n "${PATIENT_SAMPLE_LIST[$list_idx]}" \
            -R "${REFERENCE_SET_NAME}" \
            "${REGISTRY_PATH}" \
            "${DATASET_NAME}" \
            "${PATIENT_LIST[$list_idx]}" \
            "${PATIENT_SAMPLE_LIST[$list_idx]}" \
            "${GROUP}"/"${GROUP}"_${i}.vcf.gz
        echo candig_repo add-variantset -r \
            -I "${GROUP}"/"${GROUP}"_${i}.vcf.gz.tbi \
            -n "${PATIENT_SAMPLE_LIST[$list_idx]}" \
            -R "${REFERENCE_SET_NAME}" \
            "${REGISTRY_PATH}" \
            "${DATASET_NAME}" \
            "${PATIENT_LIST[$list_idx]}" \
            "${PATIENT_SAMPLE_LIST[$list_idx]}" \
            "${GROUP}"/"${GROUP}"_${i}.vcf.gz
        list_idx=$((list_idx+1))
    done