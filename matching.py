
# coding: utf-8

# In[121]:


import numpy as np
import pandas as pd
from os.path import join as opj


# ### Matching based on Volumes
# * Volume bins
#     * 100 - 150
#     * 150 - 200
#     * 200 - 250
#     * 250 - 300

# In[122]:


# ## Create a function to do volumes matching

# In[147]:


# def volumes_matching(volumes_bins, df_demographics, df_TD_phenotype, df_AUT_phenotype):
#     # Load demographics file
#
# #     demographics_file_path = '/home1/varunk/Autism-Connectome-Analysis-brain_connectivity/notebooks/demographics.csv'
# #     phenotype_file_path = '/home1/varunk/data/ABIDE1/RawDataBIDs/composite_phenotypic_file.csv'
# #     volumes_bins = np.array([[0,150],[151,200],[201,250],[251,300]])
#
#
# #     df_demographics = pd.read_csv(demographics_file_path)
#     df_demographics_volumes = df_demographics.as_matrix(['SITE_NAME','VOLUMES']).squeeze()
#
#
#
# #     df_phenotype = pd.read_csv(phenotype_file_path)
# #     df_phenotype = df_phenotype.sort_values(['SUB_ID'])
#
#
#
#     bins_volumes_AUT_data = []
#     bins_volumes_TD_data = []
#
#     for counter, _bin in enumerate(volumes_bins):
#         df_demographics_volumes_selected_bin = df_demographics_volumes[np.where(np.logical_and((df_demographics_volumes[:,1] >= _bin[0]),(df_demographics_volumes[:,1] <= _bin[1])))]
#
#
#         selected_AUT = pd.DataFrame()
#         selected_TD = pd.DataFrame()
#         for site in df_demographics_volumes_selected_bin:
# #             print(site[0])
#             selected_AUT = pd.concat([selected_AUT,df_AUT_phenotype.loc[(df_AUT_phenotype['SEX'] == 1)
#                                                                     & (df_AUT_phenotype['DSM_IV_TR'] == 1)
#                                                                     & (df_AUT_phenotype['SITE_ID'] == site[0])]])
#             selected_TD = pd.concat([selected_TD,df_TD_phenotype.loc[(df_TD_phenotype['SEX'] == 1)
#                                                                       & (df_TD_phenotype['DSM_IV_TR'] == 0)
#                                                                       & (df_TD_phenotype['SITE_ID'] == site[0])]])
#
#         bins_volumes_AUT_data.append(selected_AUT)
#         bins_volumes_TD_data.append(selected_TD)
#
#     matched_df_TD,matched_df_AUT = matching(volumes_bins, bins_volumes_TD_data, bins_volumes_AUT_data)
# #     sub_ids = selected_df_TD.as_matrix(['SUB_ID']).squeeze()
#     # matched_df_TD.to_csv('volumes_matched_TD.csv')
#     return matched_df_TD,matched_df_AUT



def matching(bins, bins_TD_data, bins_AUT_data, randomize = False):
    # num_bins = 4
    print('Original data stats')
    print('Range    ','TD ','AUT ','Ratio TD/AUT')
    ratio = np.zeros((len(bins_TD_data)))
    for i in range(len(bins_TD_data)):
        ratio[i] = bins_TD_data[i].shape[0]/bins_AUT_data[i].shape[0]
        print(bins[i],bins_TD_data[i].shape[0],bins_AUT_data[i].shape[0], ratio[i])

    min_ratio = np.min(ratio)
    min_index = np.argmin(ratio)

    new_TD = np.zeros((len(bins_TD_data)))
    new_AUT = np.zeros((len(bins_AUT_data)))

#     matched_df_AUT = None
#     matched_df_TD = None
    if min_ratio < 1:
        _ratio = 1.0 / ratio
        min_ratio = np.min(_ratio)
        print('Ratio = ',min_ratio)
#    -------------------------------------------
        print('Matched data stats')
        print('Range    ','TD ','AUT ')
        for i in range(len(bins_TD_data)):
            new_AUT[i] = np.floor(bins_TD_data[i].shape[0] * min_ratio)
            print(bins[i],bins_TD_data[i].shape[0],new_AUT[i])

        # Now loop over all the bins created and select the specific number of subjects randomly from each TD bin

        # AUT_idx_list = []
        selected_df_AUT = pd.DataFrame()
        selected_df_TD = pd.DataFrame()


        for i in range(len(bins_AUT_data)):
            idx = np.arange(len(bins_AUT_data[i]))
            if randomize == True:
                np.random.shuffle(idx)
            idx = idx[0:int(new_AUT[i])]
            # AUT_idx_list.append(idx)
            selected_df_AUT = pd.concat([selected_df_AUT, bins_AUT_data[i].iloc[idx]])
            selected_df_TD = pd.concat([selected_df_TD, bins_TD_data[i]])

        matched_df_AUT = selected_df_AUT.sort_values(['SUB_ID'])
        matched_df_TD = selected_df_TD.sort_values(['SUB_ID'])

        return matched_df_TD, matched_df_AUT


#     -------------------------------------

    print('Matched data stats')
    print('Range    ','TD ','AUT ')
    for i in range(len(bins_TD_data)):
        new_TD[i] = np.floor(bins_AUT_data[i].shape[0] * min_ratio)
        print(bins[i],new_TD[i],bins_AUT_data[i].shape[0])

    # Now loop over all the bins created and select the specific number of subjects randomly from each TD bin

    # TD_idx_list = []
    selected_df_TD = pd.DataFrame()
    selected_df_AUT = pd.DataFrame()


    for i in range(len(bins_TD_data)):
        idx = np.arange(len(bins_TD_data[i]))
        if randomize == True:
            np.random.shuffle(idx)
        idx = idx[0:int(new_TD[i])]
        # TD_idx_list.append(idx)
        selected_df_TD = pd.concat([selected_df_TD, bins_TD_data[i].iloc[idx]])
        selected_df_AUT = pd.concat([selected_df_AUT, bins_AUT_data[i]])

    matched_df_TD = selected_df_TD.sort_values(['SUB_ID'])
    matched_df_AUT = selected_df_AUT.sort_values(['SUB_ID'])

    return matched_df_TD,matched_df_AUT


# In[150]:

# Usage
# demographics_file_path = '/home1/varunk/Autism-Connectome-Analysis-brain_connectivity/notebooks/demographics.csv'
# phenotype_file_path = '/home1/varunk/data/ABIDE1/RawDataBIDs/composite_phenotypic_file.csv'
#
# df_demographics = pd.read_csv(demographics_file_path)
#
#
# df_phenotype = pd.read_csv(phenotype_file_path)
# df_phenotype = df_phenotype.sort_values(['SUB_ID'])
#
# volumes_bins = np.array([[0,150],[151,200],[201,250],[251,300]])
#
# matched_df_TD,matched_df_AUT = volumes_matching(volumes_bins, df_demographics, df_phenotype, df_phenotype)


# In[151]:


def age_matching(age_bins, df_TD_phenotype, df_AUT_phenotype, base_directory ):
#     age_bins = np.array([[0,9],[9,12],[12,15],[15,18]])

    bins_age_AUT_data = []
    bins_age_TD_data = []

    # for counter, _bin in enumerate(age_bins):

    log_path = opj(base_directory,"log.txt")
    log = open(log_path, 'a')
    log.write("------------- Age Matching with the following bins -------------\n")
    log.write("Age Bins: %s \n"%age_bins)
    log.flush()


    for age in age_bins:

        selected_AUT = pd.DataFrame()
        selected_TD = pd.DataFrame()


#         print(age[0], age[1])
        # selected_AUT = pd.concat([selected_AUT,df_AUT_phenotype[(df_AUT_phenotype['SEX'] == 1)
        #                                                         & (df_AUT_phenotype['DSM_IV_TR'] == 1)
        #                                                         & (df_AUT_phenotype['AGE_AT_SCAN'] > age[0])
        #                                                         & (df_AUT_phenotype['AGE_AT_SCAN'] <= age[1]) ]])
        # selected_TD = pd.concat([selected_TD,df_TD_phenotype.loc[(df_TD_phenotype['SEX'] == 1)
        #                                                         & (df_TD_phenotype['DX_GROUP'] == 2)
        #                                                         & (df_TD_phenotype['AGE_AT_SCAN'] > age[0])
        #                                                         & (df_TD_phenotype['AGE_AT_SCAN'] <= age[1]) ]])

        selected_AUT = pd.concat([selected_AUT,df_AUT_phenotype[(df_AUT_phenotype['AGE_AT_SCAN'] > age[0])
                                                                & (df_AUT_phenotype['AGE_AT_SCAN'] <= age[1]) ]])
        selected_TD = pd.concat([selected_TD,df_TD_phenotype.loc[(df_TD_phenotype['AGE_AT_SCAN'] > age[0])
                                                                & (df_TD_phenotype['AGE_AT_SCAN'] <= age[1]) ]])


        bins_age_AUT_data.append(selected_AUT)
        bins_age_TD_data.append(selected_TD)

    matched_df_TD,matched_df_AUT = matching(age_bins, bins_age_TD_data, bins_age_AUT_data)
#     sub_ids = selected_df_TD.as_matrix(['SUB_ID']).squeeze()
    # matched_df_TD.to_csv('age_matched_TD.csv')
    return matched_df_TD,matched_df_AUT


# In[152]:

# Usage
# demographics_file_path = '/home1/varunk/Autism-Connectome-Analysis-brain_connectivity/notebooks/demographics.csv'
# phenotype_file_path = '/home1/varunk/data/ABIDE1/RawDataBIDs/composite_phenotypic_file.csv'
#
# df_demographics = pd.read_csv(demographics_file_path)
#
#
# df_phenotype = pd.read_csv(phenotype_file_path)
# df_phenotype = df_phenotype.sort_values(['SUB_ID'])
#
# age_bins = np.array([[0,9],[9,12],[12,15],[15,18]])
#
# matched_df_TD,matched_df_AUT = age_matching(age_bins, matched_df_TD, df_phenotype)


# ## TR Matching

# In[153]:


def tr_matching(TR_bins, df_demographics, df_TD_phenotype, df_AUT_phenotype, base_directory ):
    # df_demographics = pd.read_csv(demographics_file_path)
    df_demographics_TR = df_demographics.as_matrix(['SITE_NAME','TR']).squeeze()

    log_path = opj(base_directory,"log.txt")
    log = open(log_path, 'a')
    log.write("------------- TR Matching with the following bins -------------\n")
    log.write("TR Bins: %s \n"%TR_bins)
    log.flush()



    # df_phenotype = pd.read_csv(phenotype_file_path)
    # df_phenotype = df_phenotype.sort_values(['SUB_ID'])



    bins_TR_AUT_data = []
    bins_TR_TD_data = []

    for counter, _bin in enumerate(TR_bins):
        df_demographics_TR_selected_bin = df_demographics_TR[np.where(np.logical_and((df_demographics_TR[:,1] > _bin[0]),(df_demographics_TR[:,1] <= _bin[1])))]


        selected_AUT = pd.DataFrame()
        selected_TD = pd.DataFrame()
        for site in df_demographics_TR_selected_bin:
    #             print(site[0])
            # selected_AUT = pd.concat([selected_AUT,df_AUT_phenotype.loc[(df_AUT_phenotype['SEX'] == 1)
            #                                                         & (df_AUT_phenotype['DSM_IV_TR'] == 1)
            #                                                         & (df_AUT_phenotype['SITE_ID'] == site[0])]])
            # selected_TD = pd.concat([selected_TD,df_TD_phenotype.loc[(df_TD_phenotype['SEX'] == 1)
            #                                                           & (df_TD_phenotype['DX_GROUP'] == 2)
            #                                                           & (df_TD_phenotype['SITE_ID'] == site[0])]])

            selected_AUT = pd.concat([selected_AUT,df_AUT_phenotype.loc[(df_AUT_phenotype['SITE_ID'] == site[0])]])
            selected_TD = pd.concat([selected_TD,df_TD_phenotype.loc[(df_TD_phenotype['SITE_ID'] == site[0])]])

        bins_TR_AUT_data.append(selected_AUT)
        bins_TR_TD_data.append(selected_TD)

    matched_df_TD, matched_df_AUT = matching(TR_bins, bins_TR_TD_data, bins_TR_AUT_data)
    #     sub_ids = selected_df_TD.as_matrix(['SUB_ID']).squeeze()
    # matched_df_TD.to_csv('TR_matched_TD.csv')
    return matched_df_TD, matched_df_AUT


# In[154]:


# usage
# demographics_file_path = '/home1/varunk/Autism-Connectome-Analysis-brain_connectivity/notebooks/demographics.csv'
# phenotype_file_path = '/home1/varunk/data/ABIDE1/RawDataBIDs/composite_phenotypic_file.csv'
# TR_bins = np.array([[0,2],[2,2.5],[2.5,3.0]])
#
#
# df_demographics = pd.read_csv(demographics_file_path)
# df_phenotype = pd.read_csv(phenotype_file_path)
# df_phenotype = df_phenotype.sort_values(['SUB_ID'])
#
# matched_df_TD = tr_matching(TR_bins,df_demographics, matched_df_TD, df_phenotype)


# In[155]:


# Combined Matching Usage

# demographics_file_path = '/home1/varunk/Autism-Connectome-Analysis-brain_connectivity/notebooks/demographics.csv'
# phenotype_file_path = '/home1/varunk/data/ABIDE1/RawDataBIDs/composite_phenotypic_file.csv'
# df_demographics = pd.read_csv(demographics_file_path)
# df_phenotype = pd.read_csv(phenotype_file_path)
# df_phenotype = df_phenotype.sort_values(['SUB_ID'])
#
#
#
# # Volume matching
# print('Volume Matching')
# volumes_bins = np.array([[0,150],[151,200],[201,250],[251,300]])
# matched_df_TD = df_phenotype
# matched_df_AUT = df_phenotype
# matched_df_TD, matched_df_AUT = volumes_matching(volumes_bins, df_demographics, matched_df_TD, matched_df_AUT)
#
# # TR matching
# print('TR Matching')
# TR_bins = np.array([[0,2],[2,2.5],[2.5,3.0]])
# # matched_df_TD = df_phenotype
# # matched_df_AUT = df_phenotype
# matched_df_TD,matched_df_AUT = tr_matching(TR_bins,df_demographics, matched_df_TD, matched_df_AUT)
#
#
# # Age Matching
# print('Age Matching')
# age_bins = np.array([[0,9],[9,12],[12,15],[15,18]])
# # matched_df_TD = df_phenotype
# # matched_df_AUT = df_phenotype
# matched_df_TD,matched_df_AUT = age_matching(age_bins, matched_df_TD, matched_df_AUT)
#
#
#
# matched_df_TD.loc[(matched_df_TD['SEX'] == 1) & (matched_df_TD['DSM_IV_TR'] == 0) & (matched_df_TD['EYE_STATUS_AT_SCAN'] == 2) ]
#
# matched_df_AUT.loc[(matched_df_AUT['SEX'] == 1) & (matched_df_AUT['DSM_IV_TR'] == 1) & (matched_df_AUT['EYE_STATUS_AT_SCAN'] == 2) ]
