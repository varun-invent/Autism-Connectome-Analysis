# This script select the subjects and returns the datframe of the selected participants based on certain criterias:

import pandas as pd
import numpy as np
from os.path import join as opj

def extract(phenotype_file_path, base_directory, criteria_dict):
    df_phenotype = pd.read_csv(phenotype_file_path, encoding = "ISO-8859-1")
    # df_phenotype = df_phenotype.sort_values(['SUB_ID'])

    # df_healthy = df_phenotype.loc[(df_phenotype['SEX'] == 1) & (df_phenotype['DX_GROUP'] == 2) \
    #                                                     & (df_phenotype['EYE_STATUS_AT_SCAN'] == 1) ]
    #
    # df_diseased = df_phenotype.loc[(df_phenotype['SEX'] == 1) & (df_phenotype['DSM_IV_TR'] == 1) \
    #                                                     & (df_phenotype['EYE_STATUS_AT_SCAN'] == 1) ]


    log_path = opj(base_directory,"log.txt")
    log = open(log_path, 'a')
    log.write("------------- Extracted the subjects using the following criterias -------------\n")

    for criteria, value in criteria_dict.items():
        df = df_phenotype.loc[(df_phenotype[criteria] == value)]
        df_phenotype = df
        log.write('%s == %s \n'%(criteria, value))

    log.flush()


    return df_phenotype

def extract_with_manual_query(phenotype_file_path, base_directory):
    '''
    You can change this function's query and call this instead of the extract()
    if you want to create a more complex query manually.  
    '''
    df_phenotype = pd.read_csv(phenotype_file_path)
    df_phenotype = df_phenotype.sort_values(['SUB_ID'])

    df_healthy = df_phenotype.loc[(df_phenotype['SEX'] == 1) & (df_phenotype['DX_GROUP'] == 2) \
                                                        & (df_phenotype['EYE_STATUS_AT_SCAN'] == 1) ]

    df_diseased = df_phenotype.loc[(df_phenotype['SEX'] == 1) & (df_phenotype['DSM_IV_TR'] == 1) \
                                                        & (df_phenotype['EYE_STATUS_AT_SCAN'] == 1) ]


    log_path = opj(base_directory,"log.txt")
    log = open(log_path, 'a')
    log.write("------------- Extracted the subjects using the following criterias -------------\n")

    log.write("df_healthy = df_phenotype.loc[(df_phenotype['SEX'] == 1) & (df_phenotype['DX_GROUP'] == 2) \
                                                        & (df_phenotype['EYE_STATUS_AT_SCAN'] == 1) ]\n")

    log.flush()


    return df_healthy , df_diseased


if __name__ == "__main__":
    phenotype_file_path = '/home1/varunk/data/ABIDE1/RawDataBIDs/composite_phenotypic_file.csv'
    base_directory = ''
    criteria_dict = {'SEX' : 1,'DX_GROUP' : 2, 'EYE_STATUS_AT_SCAN' : 1}
    df_healthy = extract(phenotype_file_path, base_directory, criteria_dict)
    print(df_healthy)
    df_healthy.to_csv('healthy.csv')

    df_healthy, hf_diseased = extract_with_manual_query(phenotype_file_path, base_directory)
    df.to_csv('healthy2.csv')
