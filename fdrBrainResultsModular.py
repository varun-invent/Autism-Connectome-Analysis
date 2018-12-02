
# coding: utf-8

# # Experimenting with how to do fdr correction on masked array

# In[1]:


import os
from os.path import join as opj
# from nipype.interfaces import afni
import nibabel as nib
import json
import numpy as np
import os
from os.path import join as opj

import itertools
import nibabel as nib
from multiprocessing import Pool


# In[2]:


# # Paths
#
# path_cwd = os.getcwd()
# path_split_list = path_cwd.split('/')
# s = path_split_list[0:-1] # for getting to the parent dir of pwd
# s = opj('/',*s) # *s converts list to path, # very important to add '/' in the begining so it is read as directory later
#
#
#
# # In[3]:
#
#
# # os.chdir('/home1/varunk/Autism-Connectome-Analysis-bids-related/')
# # json_path = opj(data_directory,'task-rest_bold.json')
#
# json_path = 'scripts/json/paths.json'
# with open(json_path, 'rt') as fp:
#     task_info = json.load(fp)
#
#
#
# # In[4]:
#
#
#
#
# base_directory = opj(s,task_info["base_directory_for_results"])
# motion_correction_bet_directory = task_info["motion_correction_bet_directory"]
# parent_wf_directory = task_info["parent_wf_directory"]
# # functional_connectivity_directory = task_info["functional_connectivity_directory"]
# functional_connectivity_directory = 'temp_fc'
# coreg_reg_directory = task_info["coreg_reg_directory"]
# atlas_resize_reg_directory = task_info["atlas_resize_reg_directory"]
# data_directory = opj(s,task_info["data_directory"])
# datasink_name = task_info["datasink_name"]
# # fc_datasink_name = task_info["fc_datasink_name"]
# fc_datasink_name = 'temp_dataSink'
# atlasPath = opj(s,task_info["atlas_path"])
#
# hypothesis_test_dir = opj(base_directory, task_info["hypothesis_test_dir"])





# -----------------------------------------------------------------------------------------------------------------------




# In[5]:


# base_directory


# In[6]:


# brain_voxel_list_rand = np.random.rand(10)


# In[7]:


def count_voxel_stats(pvals_list, qvals_list, map_logp_list, map_logq_list):


#     P_brain_voxel_list, Q_brain_voxel_list = Pval_Qval_tuple

    map_logp_list = np.absolute(map_logp_list)
    map_logq_list = np.absolute(map_logq_list)

#     min p value
    min_pval = np.min(pvals_list)

#     min q value
    min_qval = np.min(qvals_list)

#     p value less than 0.1
    p_lt_point_1 = np.shape(np.where(pvals_list < 0.1))[1]

#     p value less than 0.01
    p_lt_point_01 = np.shape(np.where(pvals_list < 0.01))[1]

#     p value less than 0.05
    p_lt_point_05 = np.shape(np.where(pvals_list < 0.05))[1]

#     p value less than 0.1
    q_lt_point_1 = np.shape(np.where(qvals_list < 0.1))[1]

#     p value less than 0.01
    q_lt_point_01 = np.shape(np.where(qvals_list < 0.01))[1]

#     p value less than 0.05
    q_lt_point_05 = np.shape(np.where(qvals_list < 0.05))[1]

# Voxels with abs(sign(C1MinusC2)(-1*log10(Q)))) >1.3 (t 0.5)
    logq_gt_1point3 = np.shape(np.where(map_logq_list > 1.3))[1]

# Voxels with abs(sign(C1MinusC2)(-1*log10(Q)))) >1 (t 0.1)
    logq_gt_1 = np.shape(np.where(map_logq_list > 1))[1]

# Voxels with abs(sign(C1MinusC2)(-1*log10(Q)))) >2 (t 0.01)
    logq_gt_2 = np.shape(np.where(map_logq_list > 2))[1]

# Voxels with abs(sign(C1MinusC2)(-1*log10(P)))) >1.3 (t 0.5)
    logp_gt_1point3 = np.shape(np.where(map_logp_list > 1.3))[1]

# Voxels with abs(sign(C1MinusC2)(-1*log10(P)))) >1 (t 0.1)
    logp_gt_1 = np.shape(np.where(map_logp_list > 1))[1]

# Voxels with abs(sign(C1MinusC2)(-1*log10(P)))) >2 (t 0.01)
    logp_gt_2 = np.shape(np.where(map_logp_list > 2))[1]


    return min_pval,min_qval,p_lt_point_1,p_lt_point_01,p_lt_point_05,q_lt_point_1,    q_lt_point_01,q_lt_point_05, logq_gt_1point3, logq_gt_1 ,logq_gt_2 ,logp_gt_1point3, logp_gt_1, logp_gt_2




# In[8]:


def fdr_correction_and_viz(Pvals_path, Tvals_path, C1_path, C2_path, mask_path, save_destination, affine, header, combination):
    alpha = 0.05

    Pvals = np.load(Pvals_path)
    Tvals= np.load(Tvals_path)
    C1 = np.load(C1_path)
    C2 = np.load(C2_path)



    mask = nib.load(mask_path).get_data()

    brain_indices = np.where(mask != 0 )

    from statsmodels.sandbox.stats.multicomp import fdrcorrection0

    Pvals_shape = Pvals.shape

    Qvals = np.zeros(Pvals_shape)


    map_C1MinusC2 = C1 - C2

    # sign(c1-c2) * -1 * log10(p)
    map_logp = np.multiply(np.sign(map_C1MinusC2),(-1*np.log10(Pvals)))


    roi_voxel_stats_matrix = np.zeros((Pvals_shape[3], 14)) # cozthere are 14 statistical attributes



    for roi in range(Pvals_shape[3]):

        print('Computing Stats for ROI: ',roi)

        #         pvals = ma.masked_array(Pvals[0], mask = mask, fill_value = 0)

        pvals = Pvals[:,:,:,roi]
        pvals_shape = pvals.shape

        #         inp = pvals[~pvals.mask]

        # Flatten inp and check if you get back the original matrix after
        #         inp = inp.ravel()

        pvals_list = pvals[brain_indices]

        _, qvals_list  = fdrcorrection0(pvals_list,alpha)

#       from IPython.core.debugger import Tracer; Tracer()()
        # map_logq_list = map_logq[brain_indices]
        map_logp_list = map_logp[:,:,:,roi][brain_indices]

        # print("Size of map_logp_list ",map_logp_list.shape)
#         print("Brain Indices: ", brain_indices)

        map_C1MinusC2_list = map_C1MinusC2[:,:,:,roi][brain_indices]


        #     Calculate voxel stats using the below function

        Qvals[:,:,:,roi][brain_indices] = qvals_list




        map_logq_list = np.multiply(np.sign(map_C1MinusC2_list),(-1*np.log10(qvals_list)))

        # print("Size of map_logq_list ",map_logq_list.shape)

        roi_voxel_stats_matrix[roi,:] = count_voxel_stats(pvals_list, qvals_list,                                                           map_logp_list, map_logq_list)


        # print('Stats Computed for ROI: ',roi)

    #       Save the CSV file and the Additional Brain file to visualize


    # sign(c1-c2) * -1 * log10(q)
    map_logq = np.multiply(np.sign(map_C1MinusC2),(-1*np.log10(Qvals)))



    save_destination_new = opj(save_destination,combination)
    if not os.path.exists(save_destination_new):
        os.mkdir(save_destination_new)

    print('Saving Files in directory: ', save_destination_new)

    print('Saving Stats CSV : ',)
    csv_name = 'roi_voxel_stats_' + combination + '.csv'
    np.savetxt(csv_name,roi_voxel_stats_matrix,delimiter=',',header='min_pval,min_qval,p_lt_point_1,p_lt_point_01, p_lt_point_05, q_lt_point_1, q_lt_point_01,q_lt_point_05, logq_gt_1point3, logq_gt_1 ,logq_gt_2 ,logp_gt_1point3, logp_gt_1, logp_gt_2'
              )



    print('Saving Pvals.nii.gz')
    Pvals_name = opj(save_destination_new,'Pvals.nii.gz')
    Pvals_brain_with_header = nib.Nifti1Image(Pvals, affine= affine,header = header)
    nib.save(Pvals_brain_with_header,Pvals_name)

    print('Saving Tvals.nii.gz')
    Tvals_name = opj(save_destination_new,'Tvals.nii.gz')
    Tvals_brain_with_header = nib.Nifti1Image(Tvals, affine= affine,header = header)
    nib.save(Tvals_brain_with_header,Tvals_name)

    print('Saving Qvals.nii.gz')
    Qvals_name = opj(save_destination_new,'Qvals.nii.gz')
    Qvals_brain_with_header = nib.Nifti1Image(Qvals, affine= affine,header = header)
    nib.save(Qvals_brain_with_header,Qvals_name)

    print('Saving C1MinusC2.nii.gz')
    C1MinusC2_name = opj(save_destination_new,'C1MinusC2.nii.gz')
    C1MinusC2_brain_with_header = nib.Nifti1Image(map_C1MinusC2, affine= affine,header = header)
    nib.save(C1MinusC2_brain_with_header,C1MinusC2_name)

    print('Saving map_logp.nii.gz')
    map_logp_name = opj(save_destination_new,'map_logp.nii.gz')
    map_logp_brain_with_header = nib.Nifti1Image(map_logp, affine= affine,header = header)
    nib.save(map_logp_brain_with_header,map_logp_name)

    print('Saving map_logq.nii.gz')
    map_logq_name = opj(save_destination_new,'map_logq.nii.gz')
    map_logq_brain_with_header = nib.Nifti1Image(map_logq, affine= affine,header = header)
    nib.save(map_logq_brain_with_header,map_logq_name)





# In[ ]:


# base_directory


# In[ ]:



def _main(params):
    combination, base_directory, hypothesis_test_dir, header, mask_path, affine, fdr_results_dir = params
    #
    # motion_param_regression, band_pass_filtering, global_signal_regression,
    #  smoothing, base_directory, hypothesis_test_dir, header, mask_path, affine,fdr_results_dir = params


    # combination = 'motionRegress' + str(int(motion_param_regression)) +\
    #               'global' + str(int(global_signal_regression)) + \
    #               'smoothing' + str(int(smoothing)) +\
    #               'filt' + str(int(band_pass_filtering))

    if fdr_results_dir == None:
        fdr_results_dir = 'fdr_and_results_modular'
    save_destination = opj(base_directory,fdr_results_dir,combination)
    if not os.path.exists(save_destination):
        os.makedirs(save_destination)
    os.chdir(save_destination)
    print('Saving the results in ',save_destination)

# ----------

    # motion_param_regression, band_pass_filtering, global_signal_regression,smoothing = params
    # combination = 'motionRegress' + str(int(motion_param_regression)) + 'filt' + \
    #       str(int(band_pass_filtering)) + 'global' + str(int(global_signal_regression)) + \
    #       'smoothing' + str(int(smoothing))
    print("Combination: ",combination)
    # print(motion_param_regression, band_pass_filtering, global_signal_regression,smoothing)
    Pvals_path = opj(hypothesis_test_dir,combination,'Pvals.npy')
    Tvals_path = opj(hypothesis_test_dir,combination,'Tvals.npy')
    C1_path = opj(hypothesis_test_dir,combination,'meanC1.npy')
    C2_path = opj(hypothesis_test_dir,combination,'meanC2.npy')


    fdr_correction_and_viz(Pvals_path, Tvals_path, C1_path, C2_path, mask_path,\
     save_destination, affine, header, combination )



#     print(C2_path)


def main(paths, calc_residual, smoothing, band_pass_filtering, volCorrect,
 num_proc = 7, calc_residual_options = None):
    # json_path=paths[0]
    base_directory=paths['base_directory']
    # motion_correction_bet_directory=paths[2]
    # parent_wf_directory=paths[3]
    # functional_connectivity_directory=paths[4]
    # coreg_reg_directory=paths[5]
    # atlas_resize_reg_directory=paths[6]
    # subject_list = paths[7]
    # datasink_name=paths[8]
    fc_datasink_name=paths['fc_datasink_name']
    # atlasPath=paths[10]
    brain_path=paths['brain_path']
    # mask_path=paths[12]
    # atlas_path=paths[13]
    # tr_path=paths[14]
    # motion_params_path=paths[15]
    # func2std_mat_path=paths[16]
    # MNI3mm_path=paths[17]
    # demographics_file_path = paths[18]
    # phenotype_file_path = paths[19]
    # data_directory = paths[20]
    hypothesis_test_dir = paths['hypothesis_test_dir']
    fdr_results_dir = paths['fdr_results_dir']
    binarized_atlas_mask_path = paths['binarized_atlas_mask_path']


    # import pdb; pdb.set_trace()



    # itr = (list(itertools.product([0, 1], repeat=3)))

    # itr = [(1,0,0,1)]
    # ,(1,1,1)

    # mask_path =  mni2mmMask # MNI 3mm brain voxels
    # mask_path = opj(base_directory,parent_wf_directory,motion_correction_bet_directory,coreg_reg_directory,'resample_mni/MNI152_T1_2mm_brain_resample_mask.nii.gz')
    mask_path =  binarized_atlas_mask_path
    # '/home1/varunk/atlas/Full_brain_atlas_thr0-2mm/fullbrain_atlas_thr0-3mm_binarized.nii.gz'

    # motion_param_regression, global_signal_regression, smoothing,band_pass_filtering, volCorrect = params # just execute again-- hereee

    # combination = 'motionRegress' + str(int(motion_param_regression)) +\
    #               'global' + str(int(global_signal_regression)) + \
    #               'smoothing' + str(int(smoothing)) +\
    #               'filt' + str(int(band_pass_filtering))

    comb = ''
    for a in calc_residual_options:
        comb = comb + a

    combination = 'calc_residual' + str(int(calc_residual)) + \
    'smoothing' + str(int(smoothing)) +\
    'filt' + str(int(band_pass_filtering)) +\
    'calc_residual_options' + comb

    print("Combination: ",combination)


    fc_file_list = opj(base_directory,fc_datasink_name,combination,'fc_map_brain_file_list.npy')

    # brain_path = '/home1/varunk/results_again_again/fc_motionRegress1filt1global0/_subject_id_0050002/func2std_xform/0050002_fc_map_flirt.nii.gz'

    brain_path = fc_file_list # getting a header to be used later to save brain files
    brain_path = np.load(brain_path)[0]
    brain_data = nib.load(brain_path)
    affine=brain_data.affine
    header = brain_data.header


    pool = Pool(num_proc)
    #
    # itr = [(motion_param_regression, band_pass_filtering,global_signal_regression,
    #  smoothing,base_directory, hypothesis_test_dir, header, mask_path, affine,fdr_results_dir)]

    itr = [(combination, base_directory, hypothesis_test_dir, header, mask_path, affine, fdr_results_dir)]
    data_outputs = pool.map(_main, itr)



# ------------------------------------------------------------------------------------------------------------------------
