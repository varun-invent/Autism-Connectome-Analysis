import nibabel as nib
import numpy as np
from scipy import stats
from numpy import ma
import scipy.special as special
from statsmodels.stats import multitest
import itertools
import os
from os.path import join as opj
# from nipype.interfaces import afni
import nibabel as nib
import json
import numpy as np
import matching
import pandas as pd
import ttest as tt


# Now construct a function that takes a list of SUB_ID's and returns the FC Maps
def get_subject_fc_file(subject_id_list,fc_file_path, bugs):
    import re

    return_fc_maps = []
    fc_file_list = np.load(fc_file_path)
    print('Brain files: ',fc_file_list)
    for subject_id in subject_id_list:
#         print("For subject: ",subject_id)
        found =  False
        for brain in fc_file_list:
            sub_id_extracted = re.search('.+_subject_id_(\d+)', brain).group(1)
            if str(subject_id) in bugs:
                # print("In Bugs with subject id ",subject_id)
                found = True
            elif (subject_id == int(sub_id_extracted)):
                found = True
                return_fc_maps.append(brain)
#                 print("Found for subject: ",subject_id)
        if found == False: # Some subject was not found Problem!
            print ('Unable to locate Subject: ',int(subject_id),'extracted: ',int(sub_id_extracted))
            # return 0
    return return_fc_maps




def main(paths, bugs, applyFisher, categoryInfo= None, match=1, motion_param_regression=0, global_signal_regression=0, band_pass_filtering=0, \
    smoothing=0, num_proc = 7):
    json_path=paths[0]
    base_directory=paths[1]
    motion_correction_bet_directory=paths[2]
    parent_wf_directory=paths[3]
    functional_connectivity_directory=paths[4]
    coreg_reg_directory=paths[5]
    atlas_resize_reg_directory=paths[6]
    subject_list = paths[7]
    datasink_name=paths[8]
    fc_datasink_name=paths[9]
    atlasPath=paths[10]
    brain_path=paths[11]
    mask_path=paths[12]
    atlas_path=paths[13]
    tr_path=paths[14]
    motion_params_path=paths[15]
    func2std_mat_path=paths[16]
    MNI3mm_path=paths[17]
    demographics_file_path = paths[18]
    phenotype_file_path = paths[19]
    data_directory = paths[20]
    hypothesis_test_dir = paths[21]

    #  Runall:


    if categoryInfo == None:
            # phenotype_file_path = '/home1/varunk/data/ABIDE1/RawDataBIDs/composite_phenotypic_file.csv'

        df = pd.read_csv(phenotype_file_path) # , index_col='SUB_ID'

        df = df.sort_values(['SUB_ID'])
        # df = df.sort_values(['SUB+AF8-ID'])

        if bugs == None:
            bugs = ['51232','51233','51242','51243','51244','51245','51246','51247','51270','51310','50045', '51276', '50746', '50727', '51276']

        # Bugs:
        # 50045 - ROI Missing
        # 51276, 50746, 50727 - Many in between ROIs Missing
        # 51276 - Many in between ROIs Missing

        # '0051242' in bugs


        # selecting Autistic males(DSM IV) of age <= 18 years
        # df_aut_lt18_m = df.loc[(df['SEX'] == 1) & (df['AGE_AT_SCAN'] <=18) & (df['DSM_IV_TR'] == 1) ]
        # df_aut_lt18_m = df.loc[(df['SEX'] == 1) & (df['AGE_AT_SCAN'] <=18) & (df['DSM_IV_TR'] == 1) & (df['EYE_STATUS_AT_SCAN'] == 1)] # eyes open
        # df_aut_lt18_m = df.loc[(df['SEX'] == 1) & (df['AGE_AT_SCAN'] <=18) & (df['DSM_IV_TR'] == 1) & (df['EYE_STATUS_AT_SCAN'] == 2)] # eyes closed
        # df_aut_lt18_m = df.loc[(df['SEX'] == 1) & (df['AGE_AT_SCAN'] >=12) & (df['AGE_AT_SCAN'] <=18) & (df['DSM_IV_TR'] == 1) & (df['EYE_STATUS_AT_SCAN'] == 1)] # eyes open age 12-18
        # df_aut_lt18_m = df.loc[(df['SEX'] == 1) & (df['AGE_AT_SCAN'] >=6) & (df['AGE_AT_SCAN'] <12) & (df['DSM_IV_TR'] == 1) & (df['EYE_STATUS_AT_SCAN'] == 1)] # eyes open age 6 - lt 12
        # df_aut_lt18_m = df.loc[(df['SEX'] == 1) & (df['AGE_AT_SCAN'] <=18) & (df['DSM_IV_TR'] == 1)] # AGE <= 18


        # In[205]:


        # df_aut_lt18_m.shape


        # In[206]:


        # df_td_lt18_m = df.loc[(df['SEX'] == 1) & (df['AGE_AT_SCAN'] <=18) & (df['DSM_IV_TR'] == 0) ]
        # df_td_lt18_m = df.loc[(df['SEX'] == 1) & (df['AGE_AT_SCAN'] <=18) & (df['DSM_IV_TR'] == 0) & (df['EYE_STATUS_AT_SCAN'] == 1)] # eyes open
        # df_td_lt18_m = df.loc[(df['SEX'] == 1) & (df['AGE_AT_SCAN'] <=18) & (df['DSM_IV_TR'] == 0) & (df['EYE_STATUS_AT_SCAN'] == 2)] # eyes closed
        # df_td_lt18_m = df.loc[(df['SEX'] == 1) & (df['AGE_AT_SCAN'] >=12) & (df['AGE_AT_SCAN'] <=18) & (df['DSM_IV_TR'] == 0) & (df['EYE_STATUS_AT_SCAN'] == 1)] # eyes open age 12- 18
        # df_td_lt18_m = df.loc[(df['SEX'] == 1) & (df['AGE_AT_SCAN'] >=6) & (df['AGE_AT_SCAN'] <12) & (df['DSM_IV_TR'] == 0) & (df['EYE_STATUS_AT_SCAN'] == 1)] # eyes open age 6 - lt 12
        # df_td_lt18_m = df.loc[(df['SEX'] == 1) & (df['AGE_AT_SCAN'] <=18) & (df['DSM_IV_TR'] == 0)] # AGE <= 18

        # df_aut_lt18_m = df.loc[(df['SEX'] == 1) & (df['AGE_AT_SCAN'] <=18) & (df['DSM_IV_TR'] == 0) & (df['EYE_STATUS_AT_SCAN'] == 2)] # TD eyes closed

        # In[207]:


        # df_td_lt18_m.shape


        # In[208]:


        # table_males_np = table_males.as_matrix(columns=['SUB_ID','DX_GROUP', 'DSM_IV_TR', 'AGE_AT_SCAN' ,'SEX' ,'EYE_STATUS_AT_SCAN'])


        # In[209]:






        if demographics_file_path == None:
            demographics_file_path = '/home1/varunk/Autism-Connectome-Analysis-brain_connectivity/notebooks/demographics.csv'

        if phenotype_file_path == None:
            phenotype_file_path = '/home1/varunk/data/ABIDE1/RawDataBIDs/composite_phenotypic_file.csv'


        df_demographics = pd.read_csv(demographics_file_path)
        df_phenotype = pd.read_csv(phenotype_file_path)
        df_phenotype = df_phenotype.sort_values(['SUB_ID'])



        # Volume matching
        # print('Volume Matching')
        # volumes_bins = np.array([[0,150],[151,200],[201,250],[251,300]])
        # matched_df_TD = df_phenotype
        # matched_df_AUT = df_phenotype
        # matched_df_TD, matched_df_AUT = matching.volumes_matching(volumes_bins, df_demographics, matched_df_TD, matched_df_AUT)
        #



        # Age 6 - 18 Autistic vs Healthy

        df_td_lt18_m = df_phenotype.loc[(df_phenotype['SEX'] == 1) & (df_phenotype['DX_GROUP'] == 2) \
                                                            & (df_phenotype['EYE_STATUS_AT_SCAN'] == 1) ]


        df_aut_lt18_m = df_phenotype.loc[(df_phenotype['SEX'] == 1) & (df_phenotype['DSM_IV_TR'] == 1) \
                                                            & (df_phenotype['EYE_STATUS_AT_SCAN'] == 1) ]


        # Age 6 - 18 Aspergers vs Healthy

        # df_td_lt18_m = df_phenotype.loc[(df_phenotype['SEX'] == 1) & (df_phenotype['DX_GROUP'] == 2) \
        #                                                     & (df_phenotype['EYE_STATUS_AT_SCAN'] == 1) ]
        #
        #
        # df_aut_lt18_m = df_phenotype.loc[(df_phenotype['SEX'] == 1) & (df_phenotype['DSM_IV_TR'] == 2) \
        #                                                     & (df_phenotype['EYE_STATUS_AT_SCAN'] == 1) ]


        # Age 6 - 18 Aspergers vs Autistic


        # df_td_lt18_m = df_phenotype.loc[(df_phenotype['SEX'] == 1) & (df_phenotype['DSM_IV_TR'] == 2) \
        #                                                     & (df_phenotype['EYE_STATUS_AT_SCAN'] == 1) ]
        #
        # df_aut_lt18_m = df_phenotype.loc[(df_phenotype['SEX'] == 1) & (df_phenotype['DSM_IV_TR'] == 1) \
        #                                                     & (df_phenotype['EYE_STATUS_AT_SCAN'] == 1) ]


        # ------------------  Only healthy ---------------------------------------

        # healthy_subjects = df_phenotype.loc[(df_phenotype['SEX'] == 1) & (df_phenotype['DX_GROUP'] == 2) \
        #                                                     & (df_phenotype['EYE_STATUS_AT_SCAN'] == 1) ]
        #
        #
        # # import pdb; pdb.set_trace()
        # healthy_subjects = healthy_subjects.reindex(np.random.permutation(healthy_subjects.index))
        #
        # df_td_lt18_m = healthy_subjects[0: int(healthy_subjects.shape[0]/2.0)]
        #
        # df_aut_lt18_m = healthy_subjects[int(healthy_subjects.shape[0]/2.0)+1 :]


        # ------------------  Only healthy again Saity Check---------------------------------------

        # healthy_subjects = df_phenotype.loc[(df_phenotype['SEX'] == 1) & (df_phenotype['DX_GROUP'] == 2) \
        #                                                     & (df_phenotype['EYE_STATUS_AT_SCAN'] == 1) ]
        #
        # df_td_lt18_m = healthy_subjects[0: int(healthy_subjects.shape[0]/2.0)]
        # df_aut_lt18_m = healthy_subjects[0: int(healthy_subjects.shape[0]/2.0)]

        # Age 12 - 18


        # df_td_lt18_m = matched_df_TD.loc[(matched_df_TD['SEX'] == 1) & (matched_df_TD['DSM_IV_TR'] == 0) \
        #                                                     & (matched_df_TD['EYE_STATUS_AT_SCAN'] == 1)
        #                                                     & (matched_df_TD['AGE_AT_SCAN'] >= 12 )
        #                                                     & (matched_df_TD['AGE_AT_SCAN'] <= 18) ]
        #
        # df_aut_lt18_m = matched_df_AUT.loc[(matched_df_AUT['SEX'] == 1) & (matched_df_AUT['DSM_IV_TR'] == 1) \
        #                                                     & (matched_df_AUT['EYE_STATUS_AT_SCAN'] == 1)
        #                                                     & (matched_df_AUT['AGE_AT_SCAN'] >= 12 )
        #                                                     & (matched_df_AUT['AGE_AT_SCAN'] <= 18) ]


        # --------------------- Matched data --------------------------------------------
        if match == 1:
            # # TR matching
            print('TR Matching with range (0,2.5]')

            df_demographics = df_demographics.drop(df_demographics.index[[7]]) # Deleting OHSU with volumes 82

            TR_bins = np.array([[0,2.5]])
            # TR_bins = np.array([[0,4]])


            # matched_df_TD = df_phenotype
            # matched_df_AUT = df_phenotype
            df_td_lt18_m, df_aut_lt18_m = matching.tr_matching(TR_bins,df_demographics, df_td_lt18_m, df_aut_lt18_m)


            # Age Matching
            print('Age Matching')
            age_bins = np.array([[0,9],[9,12],[12,15],[15,18]])
            # matched_df_TD = df_phenotype
            # matched_df_AUT = df_phenotype
            df_td_lt18_m, df_aut_lt18_m = matching.age_matching(age_bins, df_td_lt18_m, df_aut_lt18_m)

            df_aut = df_aut_lt18_m
            df_td = df_td_lt18_m

            df_aut_subid = df_aut.as_matrix(columns=['SUB_ID'])
            df_td_subid = df_td.as_matrix(columns=['SUB_ID'])







        # ----------------------------------------Checking the difference between eyes closed vs open ---------------------------------

        # df_td_lt18_m = matched_df_TD.loc[(matched_df_TD['SEX'] == 1) & (matched_df_TD['DSM_IV_TR'] == 0) \
        #                                                     & (matched_df_TD['EYE_STATUS_AT_SCAN'] == 1) \
        #                                                     & (matched_df_TD['AGE_AT_SCAN'] >= 12 ) \
        #                                                     & (matched_df_TD['AGE_AT_SCAN'] <= 18) ]
        #
        # df_aut_lt18_m = matched_df_TD.loc[(matched_df_TD['SEX'] == 1) & (matched_df_TD['DSM_IV_TR'] == 0) \
        #                                                     & (matched_df_TD['EYE_STATUS_AT_SCAN'] == 2) \
        #                                                     & (matched_df_TD['AGE_AT_SCAN'] >= 12 ) \
        #                                                     & (matched_df_TD['AGE_AT_SCAN'] <= 18) ]


        # -------------------------------------------------------------------------------------------------------------------------
        # import pdb; pdb.set_trace()
        # df_aut_subid = df_aut_lt18_m.as_matrix(columns=['SUB_ID'])
        # df_td_subid = df_td_lt18_m.as_matrix(columns=['SUB_ID'])
        #
        #
        # combination = 'motionRegress' + str(int(motion_param_regression)) + 'filt' + \
        #           str(int(band_pass_filtering)) + 'global' + str(int(global_signal_regression)) + \
        #           'smoothing' + str(int(smoothing))
        #
        # print("Combination: ",combination)
        # print(motion_param_regression,band_pass_filtering, global_signal_regression, smoothing)
        #
        # save_destination = opj(hypothesis_test_dir,combination)
        # print('Saving files in ',save_destination)
        # if not os.path.exists(save_destination):
        #     os.makedirs(save_destination) # to create a nested directory structure
        #
        # save_destination_TD = opj(hypothesis_test_dir,combination,'TD_subects.csv')
        # save_destination_AUT = opj(hypothesis_test_dir,combination,'AUT_subjects.csv')
        #
        # print("Storing the subjects' information used")
        # # df_td_lt18_m.to_csv('TD_subects.csv')
        # df_td_lt18_m.to_csv(save_destination_TD)
        # print('Saved TD_subects.csv')
        # # df_aut_lt18_m.to_csv('AUT_subjects.csv')
        # df_aut_lt18_m.to_csv(save_destination_AUT)
        # print('Saved AUT_subects.csv')


    # In[210]:




    else:
        df = pd.read_csv(categoryInfo) # , index_col='SUB_ID'
        df_phenotype = df.sort_values(['PART_ID'])

        df_aut = df_phenotype.loc[(df_phenotype['Control'] == 1)]
        df_td = df_phenotype.loc[(df_phenotype['Control'] == 0)]

        df_aut_subid = df_aut.as_matrix(columns=['PART_ID'])
        df_td_subid = df_td.as_matrix(columns=['PART_ID'])

    # import pdb; pdb.set_trace()


    combination = 'motionRegress' + str(int(motion_param_regression)) + \
     'global' + str(int(global_signal_regression)) + 'smoothing' + str(int(smoothing)) +\
     'filt' + str(int(band_pass_filtering))

    print("Combination: ",combination)
    print(motion_param_regression,band_pass_filtering, global_signal_regression, smoothing)

    save_destination = opj(hypothesis_test_dir,combination)
    print('Saving files in ',save_destination)
    if not os.path.exists(save_destination):
        os.makedirs(save_destination) # to create a nested directory structure

    save_destination_TD = opj(hypothesis_test_dir,combination,'TD_subects.csv')
    save_destination_AUT = opj(hypothesis_test_dir,combination,'AUT_subjects.csv')

    print("Storing the subjects' information used")
    # df_td_lt18_m.to_csv('TD_subects.csv')
    df_td.to_csv(save_destination_TD)
    print('Saved TD_subects.csv')
    # df_aut_lt18_m.to_csv('AUT_subjects.csv')
    df_aut.to_csv(save_destination_AUT)
    print('Saved AUT_subects.csv')




    # for motion_param_regression, band_pass_filtering, global_signal_regression, smoothing in itr:

    fc_file_list = opj(base_directory,fc_datasink_name,combination,'fc_map_brain_file_list.npy')

    print('Reading the brain paths from: ',fc_file_list)
#     apply_fisher = True

    # import pdb;pdb.set_trace()
    autistic_list = (get_subject_fc_file(df_aut_subid.squeeze(), fc_file_list, bugs))
    print("Number of autistic participants ", len(autistic_list))

    td_list = (get_subject_fc_file(df_td_subid.squeeze(), fc_file_list, bugs))
    print("Number of TD participants ", len(td_list))

    # participants_considered = min(len(autistic_list), len(td_list))

    # participants_considered = 2

    # print("Number of participants being Considered per group:", participants_considered)

    autistic_list = autistic_list#[0:participants_considered]
    td_list = td_list#[0:participants_considered]

    # Created the below mask manually using BET
    # mask = opj(base_directory,parent_wf_directory,motion_correction_bet_directory,coreg_reg_directory,'resample_mni/MNI152_T1_2mm_brain_resample_mask.nii.gz')
    # mask = opj(base_directory,parent_wf_directory,motion_correction_bet_directory,coreg_reg_directory,'atlas_resize_reg_directory/resample_atlas/fullbrain_atlas_thr0-2mm_resample_binarize.nii.gz')
    mask = '/home1/varunk/atlas/Full_brain_atlas_thr0-2mm/fullbrain_atlas_thr0-3mm_binarized.nii.gz'

    # print('Saving the results in ', hypothesis_test_dir)
    tt.main(autistic_list,td_list, combination, mask, applyFisher, hypothesis_test_dir) # sends the file path of autistic and TD and processing params
