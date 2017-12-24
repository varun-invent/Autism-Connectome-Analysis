
# coding: utf-8

# # Hypothesis Testing
# This code does the following:
# * Reads the FC files for all the subjects
# * Z-Standardize all the voxel-roi correlation values of each ROI
# * Perform two tailed t-test for each voxel-roi pair correlation across subjects (Autism vs TD)

# In[192]:


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


# In[193]:


# Paths

path_cwd = os.getcwd()
path_split_list = path_cwd.split('/')
s = path_split_list[0:-1] # for getting to the parent dir of pwd
s = opj('/',*s) # *s converts list to path, # very important to add '/' in the begining so it is read as directory later



# In[194]:



# json_path = opj(data_directory,'task-rest_bold.json')

json_path = 'scripts/json/paths.json'
with open(json_path, 'rt') as fp:
    task_info = json.load(fp)



# In[195]:




base_directory = opj(s,task_info["base_directory_for_results"])
motion_correction_bet_directory = task_info["motion_correction_bet_directory"]
parent_wf_directory = task_info["parent_wf_directory"]
# functional_connectivity_directory = task_info["functional_connectivity_directory"]
functional_connectivity_directory = 'temp_fc'
coreg_reg_directory = task_info["coreg_reg_directory"]
atlas_resize_reg_directory = task_info["atlas_resize_reg_directory"]
data_directory = opj(s,task_info["data_directory"])
datasink_name = task_info["datasink_name"]
fc_datasink_name = task_info["fc_datasink_name"]
# fc_datasink_name = 'temp_dataSink'
atlasPath = opj(s,task_info["atlas_path"])

hypothesis_test_dir = opj(base_directory, task_info["hypothesis_test_dir"])



# In[196]:


# brain_path = opj(base_directory,datasink_name,'preprocessed_brain_paths/brain_file_list.npy')
# mask_path = opj(base_directory,datasink_name,'preprocessed_mask_paths/mask_file_list.npy')
# atlas_path = opj(base_directory,datasink_name,'atlas_paths/atlas_file_list.npy')
# tr_path = opj(base_directory,datasink_name,'tr_paths/tr_list.npy')
# motion_params_path = opj(base_directory,datasink_name,'motion_params_paths/motion_params_file_list.npy')
# func2std_mat_path = opj(base_directory, datasink_name,'joint_xformation_matrix_paths/joint_xformation_matrix_file_list.npy')

# # brain_list = np.load('../results_again_again/ABIDE1_Preprocess_Datasink/preprocessed_brain_paths/brain_file_list.npy')


# In[197]:


# brain_path,mask_path,atlas_path,tr_path,motion_params_path,func2std_mat_path


# In[198]:


# brain_path = np.load(brain_path)
# mask_path = np.load(mask_path)
# atlas_path = np.load(atlas_path)
# tr_path = np.load(tr_path)
# motion_params_path = np.load(motion_params_path)


# In[199]:












#  Runall:
import pandas as pd

df = pd.read_csv('/home1/varunk/data/ABIDE1/RawDataBIDs/composite_phenotypic_file.csv') # , index_col='SUB_ID'

df = df.sort_values(['SUB_ID'])
# df = df.sort_values(['SUB+AF8-ID'])

bugs = ['51232','51233','51242','51243','51244','51245','51246','51247','51270','51310','50045', '51276', '50746', '50727', '51276']
# Bugs:
# 50045 - ROI Missing
# 51276, 50746, 50727 - Many in between ROIs Missing
# 51276 - Many in between ROIs Missing


# In[202]:


# '0051242' in bugs


# In[203]:


# df


# In[204]:


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


# --------------------- Matched data --------------------------------------------

demographics_file_path = '/home1/varunk/Autism-Connectome-Analysis-brain_connectivity/notebooks/demographics.csv'
phenotype_file_path = '/home1/varunk/data/ABIDE1/RawDataBIDs/composite_phenotypic_file.csv'
df_demographics = pd.read_csv(demographics_file_path)
df_phenotype = pd.read_csv(phenotype_file_path)
df_phenotype = df_phenotype.sort_values(['SUB_ID'])



# Volume matching
print('Volume Matching')
volumes_bins = np.array([[0,150],[151,200],[201,250],[251,300]])
matched_df_TD = df_phenotype
matched_df_AUT = df_phenotype
matched_df_TD, matched_df_AUT = matching.volumes_matching(volumes_bins, df_demographics, matched_df_TD, matched_df_AUT)

# TR matching
print('TR Matching')
TR_bins = np.array([[0,2],[2,2.5],[2.5,3.0]])
# matched_df_TD = df_phenotype
# matched_df_AUT = df_phenotype
matched_df_TD,matched_df_AUT = matching.tr_matching(TR_bins,df_demographics, matched_df_TD, matched_df_AUT)


# Age Matching
print('Age Matching')
age_bins = np.array([[0,9],[9,12],[12,15],[15,18]])
# matched_df_TD = df_phenotype
# matched_df_AUT = df_phenotype
matched_df_TD,matched_df_AUT = matching.age_matching(age_bins, matched_df_TD, matched_df_AUT)



df_td_lt18_m = matched_df_TD.loc[(matched_df_TD['SEX'] == 1) & (matched_df_TD['DSM_IV_TR'] == 0) \
                                                    & (matched_df_TD['EYE_STATUS_AT_SCAN'] == 1) ]

df_aut_lt18_m = matched_df_AUT.loc[(matched_df_AUT['SEX'] == 1) & (matched_df_AUT['DSM_IV_TR'] == 1) \
                                                    & (matched_df_AUT['EYE_STATUS_AT_SCAN'] == 1) ]

# ----------------------------------------Checking the difference between eyes closed vs open ---------------------------------
# 
# df_td_lt18_m = matched_df_TD.loc[(matched_df_TD['SEX'] == 1) & (matched_df_TD['DSM_IV_TR'] == 0) \
#                                                     & (matched_df_TD['EYE_STATUS_AT_SCAN'] == 1) ]
#
# df_aut_lt18_m = matched_df_AUT.loc[(matched_df_AUT['SEX'] == 1) & (matched_df_AUT['DSM_IV_TR'] == 0) \
#                                                     & (matched_df_AUT['EYE_STATUS_AT_SCAN'] == 2) ]
#


# -------------------------------------------------------------------------------------------------------------------------
# import pdb; pdb.set_trace()
df_aut_subid = df_aut_lt18_m.as_matrix(columns=['SUB_ID'])
df_td_subid = df_td_lt18_m.as_matrix(columns=['SUB_ID'])

print("Storing the subjects' information used")
df_td_lt18_m.to_csv('TD_subects.csv')
print('Saved TD_subects.csv')
df_aut_lt18_m.to_csv('AUT_subjects.csv')
print('Saved AUT_subects.csv')


# In[210]:


# df_aut_subid#, df_td_subid


# In[211]:


# Now construct a function that takes a list of SUB_ID's and returns the FC Maps
def get_subject_fc_file(subject_id_list,fc_file_path, bugs):
    import re

    return_fc_maps = []
    fc_file_list = np.load(fc_file_path)
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
            return 0
    return return_fc_maps





# In[212]:





# In[213]:





# In[214]:


# len(autistic_list),len(td_list)


# In[215]:


# # To Stop execution Raise error:
# raise Exception('Execution stops here!')


# In[216]:


# number_of_fcmaps = len(fc_file_list) #184


# In[217]:


# number_of_fcmaps


# In[218]:


# Author Deepak Singla: singlakdeepak5@gmail.com


def div0( a, b ):
    '''
    It is meant for ignoring the places where standard deviation
    is zero.
    '''
    """ ignore / 0, div0( [-1, 0, 1], 0 ) -> [0, 0, 0] """
    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.divide( a, b )
        # c[ ~ np.isfinite( c )] = 0  # -inf inf NaN
    return c

def calc_mean_and_std(ROICorrMaps, n_subjects, ROIAtlasmask, ddof =1, applyFisher = False):
    '''
	Function for calculating the mean and standard
	deviation of the data. At a time, only one of the nii
	file is loaded and the elements keep on adding as we
	enumerate over the subjects.
	'''
    mask = nib.load(ROIAtlasmask).get_data()
    mask = ma.masked_object(mask,0).mask
    if (n_subjects != 0):
        f = nib.load(ROICorrMaps[0])
        dimensions = f.get_header().get_data_shape()
        print(dimensions)
    else:
        exit
    mask = np.repeat(mask[:, :, :, np.newaxis], dimensions[3], axis=3)
#     print(ROICorrMaps)

    Sample_mean_Array = np.zeros(dimensions)
    Sample_std_Array = np.zeros(dimensions)
    Sample_mean_Array = ma.masked_array(Sample_mean_Array,
                                       mask = mask,
                                        fill_value = 0)
    Sample_std_Array = ma.masked_array(Sample_std_Array,
                                      mask = mask ,
                                       fill_value = 0)
    for count, subject in enumerate(ROICorrMaps):

        Corr_data = nib.load(subject).get_data()
        Corr_data = ma.masked_array(Corr_data, mask = mask, fill_value = 0)
        if applyFisher:
            Corr_data = np.arctanh(Corr_data)

        Sample_mean_Array += Corr_data
        Sample_std_Array += np.square(Corr_data)
        print('Done subject ', count+1)
    Sample_mean_Array /= n_subjects
    # import pdb; pdb.set_trace()
    Sample_std_Array = np.sqrt((Sample_std_Array - n_subjects*np.square(Sample_mean_Array))/(n_subjects - ddof))
    return Sample_mean_Array,Sample_std_Array

def calc_mean_and_std_if_npy(ROICorrMaps, n_subjects, ddof =1, applyFisher = False):
    '''
    Function to be used if the file is given in the format
    No of ROIs versus All brain voxels in the ROI mapped.
    '''
    print(ROICorrMaps)
    initialize = np.load(ROICorrMaps[0])
    initialize = ma.masked_array(initialize)
    if applyFisher:
        initialize = np.arctanh(initialize)
    Sample_mean_Array = ma.masked_array(initialize,
                                        fill_value = 0)
    Sample_std_Array = ma.masked_array(np.square(initialize),
                                       fill_value = 0)
    del initialize
    print('Done subject ', 0)
    for count, subject in enumerate(ROICorrMaps[1:]):

        Corr_data = np.load(subject)
        Corr_data = ma.masked_array(Corr_data)
        if applyFisher:
            Corr_data = np.arctanh(Corr_data)
        Sample_mean_Array += Corr_data
        Sample_std_Array += np.square(Corr_data)
        print('Done subject ', count+1)
    Sample_mean_Array /= n_subjects
    Sample_std_Array = np.sqrt((Sample_std_Array - n_subjects*np.square(Sample_mean_Array))/(n_subjects - ddof))
    return Sample_mean_Array,Sample_std_Array

def _ttest_1samp(Sample_mean_Array, Sample_std_Array, n_subjects, PopMean = 0.0):
    ttest_1samp_for_all = div0((Sample_mean_Array - PopMean)                             * np.sqrt(n_subjects), Sample_std_Array)
    df = n_subjects - 1
    # pval = stats.t.sf(np.abs(ttest_1samp_for_all), df)*2
    pval = special.betainc(0.5*df, 0.5, df/             (df + ttest_1samp_for_all*ttest_1samp_for_all)).reshape(ttest_1samp_for_all.shape)
    # ttest_1samp_for_all, pval = ma.filled(ttest_1samp_for_all), ma.filled(pval)

    return ttest_1samp_for_all, pval


def ttest_1samp_for_all_ROIs(ROICorrMaps,
                                ROIAtlasmask,
                                PopMean = 0.0,
                                applyFisher = False):
    '''
    This is the 1 sample t-test for ROI correlation maps.
    df = no of subjects - 1
    * ROICorrMaps is the list of filepaths of ROI correlation
    maps for a group.
    * Each ROI correlation map has the 4th dimension equal to
    the number of ROIs.
    * It calculates both the ttest as well as the p values.

    QUESTIONS???????????????????????????????????????????????
    For application of the Fisher transform, I saw that it is
    same as the inverse hyperbolic tangent function.
    Doubt is regarding the standard deviation of the distribution after
    applying Fisher. It was written that the sd is now 1/sqrt(no_of_subjs - 3).
    So, that means for each voxel or variable, the sd now becomes this.

    Ref: https://docs.scipy.org/doc/numpy/reference/generated/numpy.arctanh.html
     https://en.wikipedia.org/wiki/Fisher_transformation
    TO BE ASKED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # The ttest will return t value Inf or NaN where the denom is
    # zero. See what to return in these places. Ask tomorrow.
    '''


    n_subjects = len(ROICorrMaps)
    assert (n_subjects>0)
    Sample_mean_Array, Sample_std_Array = calc_mean_and_std(ROICorrMaps,
                                                            n_subjects,
                                                            ROIAtlasmask, ddof =1,
                                                            applyFisher = applyFisher)
    ttest_1samp_for_all, pval = _ttest_1samp(Sample_mean_Array,
                                             Sample_std_Array,
                                             n_subjects,
                                             PopMean = PopMean)
    return ttest_1samp_for_all, pval


def ttest_1samp_ROIs_if_npy(ROICorrMaps,
                            PopMean = 0.0,
                            applyFisher = False):
    n_subjects = len(ROICorrMaps)
    assert (n_subjects>0)
    Sample_mean_Array, Sample_std_Array =                                 calc_mean_and_std_if_npy( ROICorrMaps,
                                                        n_subjects, ddof =1,
                                                        applyFisher = applyFisher)
    return _ttest_1samp(Sample_mean_Array,
                        Sample_std_Array,
                        n_subjects,
                        PopMean = PopMean)


def _ttest_ind(Sample_mean_ArrayA, Sample_var_ArrayA, n_subjectsA,
                Sample_mean_ArrayB,Sample_var_ArrayB, n_subjectsB,
                equal_var = True):
    if equal_var:
        # force df to be an array for masked division not to throw a warning
        df = ma.asanyarray(n_subjectsA + n_subjectsB - 2.0)
        svar = ((n_subjectsA-1)*Sample_var_ArrayA+(n_subjectsB-1)*Sample_var_ArrayB)/ df
        denom = ma.sqrt(svar*(1.0/n_subjectsA + 1.0/n_subjectsB))  # n-D computation here!
    else:
        vn1 = Sample_var_ArrayA/n_subjectsA
        vn2 = Sample_var_ArrayB/n_subjectsB
        df = (vn1 + vn2)**2 / (vn1**2 / (n_subjectsA - 1) + vn2**2 / (n_subjectsB - 1))

        # If df is undefined, variances are zero.
        # It doesn't matter what df is as long as it is not NaN.
        df = np.where(np.isnan(df), 1, df)
        denom = ma.sqrt(vn1 + vn2)

    with np.errstate(divide='ignore', invalid='ignore'):
        ttest_ind = (Sample_mean_ArrayA - Sample_mean_ArrayB) / denom
    pvalues = special.betainc(0.5*df, 0.5, df/(df + ttest_ind*ttest_ind)).reshape(ttest_ind.shape)

    # ttest_ind, pvalues = ma.filled(ttest_ind), ma.filled(pvalues)
    return ttest_ind, pvalues,Sample_mean_ArrayA ,Sample_mean_ArrayB

def ttest_ind_samples(ROICorrMapsA, ROICorrMapsB, ROIAtlasmask,
                    equal_var = True, applyFisher = False):
    '''
    Modified from https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.stats.ttest_ind.html ,
    https://github.com/scipy/scipy/blob/v0.19.1/scipy/stats/stats.py#L3950-L4072

    Since it didn't support if the data is large and everything can't be loaded at once. So,
    such modification has been made.
    '''

    n_subjectsA = len(ROICorrMapsA)
    n_subjectsB = len(ROICorrMapsB)
    assert (n_subjectsA > 0)
    assert (n_subjectsB > 0)
    Sample_mean_ArrayA, Sample_std_ArrayA = calc_mean_and_std(ROICorrMapsA,
                                                              n_subjectsA,
                                                              ROIAtlasmask, ddof =1,
                                                              applyFisher = applyFisher)
    Sample_var_ArrayA = np.square(Sample_std_ArrayA)
    del(Sample_std_ArrayA)

    # n_subjectsB = len(ROICorrMapsB)
    Sample_mean_ArrayB, Sample_std_ArrayB = calc_mean_and_std(ROICorrMapsB,
                                                              n_subjectsB,
                                                              ROIAtlasmask, ddof =1,
                                                              applyFisher = applyFisher)
    Sample_var_ArrayB = np.square(Sample_std_ArrayB)
    del(Sample_std_ArrayB)

    # pvalues = stats.t.sf(np.abs(ttest_ind), df)*2
    return _ttest_ind(Sample_mean_ArrayA, Sample_var_ArrayA, n_subjectsA,
                Sample_mean_ArrayB, Sample_var_ArrayB, n_subjectsB,
                equal_var = equal_var)

def ttest_ind_samples_if_npy(ROICorrMapsA, ROICorrMapsB, equal_var = True, applyFisher = False):
    '''
    Modified from https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.stats.ttest_ind.html ,
    https://github.com/scipy/scipy/blob/v0.19.1/scipy/stats/stats.py#L3950-L4072

    Since it didn't support if the data is large and everything can't be loaded at once. So,
    such modification has been made.
    '''

    n_subjectsA = len(ROICorrMapsA)
    n_subjectsB = len(ROICorrMapsB)
    assert (n_subjectsA > 0)
    assert (n_subjectsB > 0)
    Sample_mean_ArrayA, Sample_std_ArrayA = calc_mean_and_std_if_npy(ROICorrMapsA,
                                                              n_subjectsA, ddof = 1,
                                                              applyFisher = applyFisher)
    Sample_var_ArrayA = np.square(Sample_std_ArrayA)
    del(Sample_std_ArrayA)
    Sample_mean_ArrayB, Sample_std_ArrayB = calc_mean_and_std_if_npy(ROICorrMapsB,
                                                              n_subjectsB, ddof =1,
                                                              applyFisher = applyFisher)
    Sample_var_ArrayB = np.square(Sample_std_ArrayB)
    del(Sample_std_ArrayB)
    # pvalues = stats.t.sf(np.abs(ttest_ind), df)*2
    return _ttest_ind(Sample_mean_ArrayA, Sample_var_ArrayA, n_subjectsA,
                Sample_mean_ArrayB, Sample_var_ArrayB, n_subjectsB,
                equal_var = equal_var)

def convert_ma_to_np(MaskedArrayObj):
    return ma.filled(MaskedArrayObj)

def fdr_correction(pvalues , type = 'ind_ROIs'):
    '''
    Two types:
    ind_ROIs: When the ROIs are taken independently and the FDR is done considering the
           the tests only in that ROI.
    all: When all the tests are treated as one.
    '''



# ### Create an MNI 3mm brain mask
#

# In[219]:


mask = opj(base_directory,parent_wf_directory,motion_correction_bet_directory,coreg_reg_directory,'resample_mni/MNI152_T1_2mm_brain_resample_mask.nii.gz')


# In[220]:

#  Author Deepak Singla : singlakdeepak5@gmail.com


def main_test(autistic_list, td_list,combination):


    combination = 'hypothesis_test_' + combination
    apply_fisher = True
    list1 = autistic_list
    list2 = td_list
    Tvals, Pvals, meanC1, meanC2 = ttest_ind_samples(list1,list2,mask,equal_var = False, applyFisher=apply_fisher)
    Tvals, Pvals, meanC1, meanC2 = convert_ma_to_np(Tvals), convert_ma_to_np(Pvals), convert_ma_to_np(meanC1), convert_ma_to_np(meanC2)
    save_destination = opj(hypothesis_test_dir,combination)
    if not os.path.exists(save_destination):
        os.makedirs(save_destination) # to create a nested directory structure
    Tvals_path = opj(save_destination,'Tvals')
    Pvals_path = opj(save_destination,'Pvals')
    mean1_path = opj(save_destination,'meanC1')
    mean2_path = opj(save_destination,'meanC2')
    np.save(Tvals_path,Tvals)
    np.save(Pvals_path,Pvals)
    np.save(mean1_path,meanC1)
    np.save(mean2_path,meanC2)



#
# fc_datasink_name = 'fc_datasink'
# itr = (list(itertools.product([0, 1], repeat=3)))
# (1,1,1),
itr = [(1,1,0,1)]

for motion_param_regression, band_pass_filtering, global_signal_regression, smoothing in itr:
    combination = 'pearcoff_motionRegress' + str(int(motion_param_regression)) + 'filt' + \
              str(int(band_pass_filtering)) + 'global' + str(int(global_signal_regression)) + \
              'smoothing' + str(int(smoothing))

    print("Combination: ",combination)
    print(motion_param_regression,band_pass_filtering, global_signal_regression, smoothing)
    fc_file_list = opj(base_directory,fc_datasink_name,combination,'fc_map_brain_file_list.npy')
#     print(fc_file_list)
#     apply_fisher = True


    autistic_list = (get_subject_fc_file(df_aut_subid.squeeze(), fc_file_list, bugs))
    print("Number of autistic participants ", len(autistic_list))

    td_list = (get_subject_fc_file(df_td_subid.squeeze(), fc_file_list, bugs))
    print("Number of TD participants ", len(td_list))

    # participants_considered = min(len(autistic_list), len(td_list))

    # participants_considered = 2

    # print("Number of participants being Considered per group:", participants_considered)

    autistic_list = autistic_list#[0:participants_considered]
    td_list = td_list#[0:participants_considered]

    main_test(autistic_list,td_list, combination)
