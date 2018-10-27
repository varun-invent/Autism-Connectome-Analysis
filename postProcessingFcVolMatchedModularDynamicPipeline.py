
# coding: utf-8

# ## Post processing:
# * Global Signal Regression using orthogonalization
# * Band Pass filtering 0.1 - 0.01 Hz
# * Motion regression using GLM
#
#

# In[460]:



from nipype.interfaces.fsl import (BET, ExtractROI, FAST, FLIRT, ImageMaths,
                                   MCFLIRT, SliceTimer, Threshold,Info, ConvertXFM,MotionOutliers)
from nipype.interfaces.afni import Resample
from nipype.interfaces.io import DataSink
from nipype.pipeline import Node, MapNode, Workflow, JoinNode
from nipype.interfaces.utility import IdentityInterface, Function
import os
from os.path import join as opj
from nipype.interfaces import afni
import nibabel as nib
import json
import numpy as np
import pandas as pd
from confounds import confounds_creation.calc_residuals as func_calc_residuals




# Paths

# path_cwd = os.getcwd()
# path_split_list = path_cwd.split('/')
# s = path_split_list[0:-1] # for getting to the parent dir of pwd
# s = opj('/',*s) # *s converts list to path, # very important to add '/' in the begining so it is read as directory later
#
#
#
# # In[462]:
#
#
#
# # json_path = opj(data_directory,'task-rest_bold.json')
#
# json_path = 'scripts/json/paths.json'
# with open(json_path, 'rt') as fp:
#     task_info = json.load(fp)



#
# base_directory = opj(s,task_info["base_directory_for_results"])
# motion_correction_bet_directory = task_info["motion_correction_bet_directory"]
# parent_wf_directory = task_info["parent_wf_directory"]
# functional_connectivity_directory = task_info["functional_connectivity_directory"]
# # functional_connectivity_directory = 'temp_fc'
# coreg_reg_directory = task_info["coreg_reg_directory"]
# atlas_resize_reg_directory = task_info["atlas_resize_reg_directory"]
# data_directory = opj(s,task_info["data_directory"])
# datasink_name = task_info["datasink_name"]
# fc_datasink_name = task_info["fc_datasink_name"]
# # fc_datasink_name = 'temp_dataSink'
# atlasPath = opj(s,task_info["atlas_path"])
#
#
#
# brain_path = opj(base_directory,datasink_name,'preprocessed_brain_paths/brain_file_list.npy')
# mask_path = opj(base_directory,datasink_name,'preprocessed_mask_paths/mask_file_list.npy')
# atlas_path = opj(base_directory,datasink_name,'atlas_paths/atlas_file_list.npy')
# tr_path = opj(base_directory,datasink_name,'tr_paths/tr_list.npy')
# motion_params_path = opj(base_directory,datasink_name,'motion_params_paths/motion_params_file_list.npy')
#
# func2std_mat_path = opj(base_directory, datasink_name,'joint_xformation_matrix_paths/joint_xformation_matrix_file_list.npy')
#
# MNI3mm_path = opj(base_directory,parent_wf_directory,motion_correction_bet_directory,coreg_reg_directory,'resample_mni/MNI152_T1_2mm_brain_resample.nii')

def main(paths, vols, motion_param_regression=0, global_signal_regression=0, band_pass_filtering=0, \
smoothing=0, volcorrect = 0, number_of_skipped_volumes=4, num_proc = 7, save_npy = 0 ):
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


    print('Brain Paths',brain_path)
    print('Atlas Path', atlas_path)


    brain_path = np.load(brain_path)
    mask_path = np.load(mask_path)
    atlas_path = np.load(atlas_path)
    tr_path = np.load(tr_path)
    motion_params_path = np.load(motion_params_path)
    func2std_mat_path = np.load(func2std_mat_path)




    # layout = BIDSLayout(data_directory)
    #
    # # number_of_subjects = 2 # Number of subjects you wish to preprocess
    #
    # if number_of_subjects == -1:
    #     number_of_subjects = len(layout.get_subjects())
    #
    #
    # subject_list = (layout.get_subjects())[0:number_of_subjects]
    subject_list = list(map(int, subject_list))

    if volcorrect == 1:
        subject_list, subid_vol_dict = volumeCorrect(subject_list, phenotype_file_path, demographics_file_path, vols )
    else:
        subid_vol_dict = None


    _main(subject_list, vols, subid_vol_dict, number_of_skipped_volumes, brain_path,\
        mask_path,\
        atlas_path,\
        tr_path,\
        motion_params_path,\
        func2std_mat_path,\
        MNI3mm_path,\
        base_directory,\
        fc_datasink_name,\
        motion_param_regression,\
        band_pass_filtering,\
        global_signal_regression,\
        smoothing,\
        volcorrect,\
        num_proc,\
        functional_connectivity_directory,save_npy )



def volumeCorrect(subject_list, phenotype_file_path = None, demographics_file_path = None, vols = None):
    if demographics_file_path == None:
        raise Exception('Demographics file not supplied')
        # demographics_file_path = '/home1/varunk/Autism-Connectome-Analysis-brain_connectivity/notebooks/demographics.csv'
    if phenotype_file_path == None:
        raise Exception('Phenotype file not supplied')
        # phenotype_file_path = '/home1/varunk/data/ABIDE1/RawDataBIDs/composite_phenotypic_file.csv'

    df_phenotype = pd.read_csv(phenotype_file_path)

    df_phenotype = df_phenotype.sort_values(['SUB_ID'])
    df_phenotype_sub_id = df_phenotype.as_matrix(['SITE_ID','SUB_ID']).squeeze()

    df_demographics = pd.read_csv(demographics_file_path)
    df_demographics_volumes = df_demographics.as_matrix(['SITE_NAME','VOLUMES']).squeeze()


    # SUB_ID - Volumes Dictionary
    site_vol_dict = dict(zip(df_demographics_volumes[:,0], df_demographics_volumes[:,1]))

    subid_vol_dict = dict(zip(df_phenotype_sub_id[:,1],[site_vol_dict[site] for site in df_phenotype_sub_id[:,0]] ))


    # Delete the site index which has volume < vols

    if vols == None:
        print('Performing volume correction with DEFAULT number of volumes: %s '%vols)
        vols = 120


    del_idx = []
    for idx,df in enumerate(df_demographics_volumes):
    #     print(idx,df[1])
        if df[1] < vols:
            del_idx.append(idx)

    df_demographics_volumes = np.delete(df_demographics_volumes,del_idx, axis = 0)



    df_demographics_sites_refined = df_demographics_volumes[:,0]



    subjects_refined = []
    for df in df_phenotype_sub_id:
        if df[0] in df_demographics_sites_refined:
    #         print(df[1])
            subjects_refined.append(df[1])

    number_of_subjects = len(subject_list)


    subjects_refined = list(set(subjects_refined) - (set(df_phenotype_sub_id[:,1]) - set(subject_list) ) )



    subject_list = subjects_refined[0:number_of_subjects]

    return subject_list, subid_vol_dict

def _main(subject_list,vols,subid_vol_dict, number_of_skipped_volumes,brain_path,\
    mask_path,\
    atlas_path,\
    tr_path,\
    motion_params_path,\
    func2std_mat_path,\
    MNI3mm_path,\
    base_directory,\
    fc_datasink_name,\
   motion_param_regression,\
   band_pass_filtering,\
   global_signal_regression,\
   smoothing,\
   volcorrect,\
   num_proc,\
   functional_connectivity_directory, save_npy ):

    # ## Volume correction
    # * I have already extracted 4 volumes.
    # * Now extract 120 - 4 = 116 volumes from each subject
    # * So define vols = 114
    #

    if number_of_skipped_volumes == None:
        number_of_skipped_volumes = 4
    vols = vols - number_of_skipped_volumes

    def vol_correct(sub_id, subid_vol_dict, vols, number_of_skipped_volumes):
        sub_vols = subid_vol_dict[sub_id] - number_of_skipped_volumes
        if sub_vols > vols:
            t_min = sub_vols - vols
        elif sub_vols == vols:
            t_min = 0
        else:
            raise Exception('Volumes of Sub ',sub_id,' less than desired!')
        return int(t_min)


    # In[491]:



    volCorrect = Node(Function(function=vol_correct, input_names=['sub_id','subid_vol_dict','vols','number_of_skipped_volumes'],
                                    output_names=['t_min']), name='volCorrect')

    volCorrect.inputs.subid_vol_dict = subid_vol_dict
    volCorrect.inputs.vols = vols
    volCorrect.inputs.number_of_skipped_volumes = number_of_skipped_volumes


    # ## Define a function to fetch the filenames of a particular subject ID



    def get_subject_filenames(subject_id,brain_path,mask_path,atlas_path,tr_path,motion_params_path,func2std_mat_path,MNI3mm_path):
        import re
        from itertools import zip_longest
        for brain,mask,atlas,tr,motion_param,func2std_mat in zip_longest(brain_path,mask_path,atlas_path,tr_path,motion_params_path,func2std_mat_path): #itertools helps to zip unequal save_file_list_in_mask
        #  Source : https://stackoverflow.com/questions/11318977/zipping-unequal-lists-in-python-in-to-a-list-which-does-not-drop-any-element-fro
            print('*******************',brain,mask,atlas,tr,motion_param,func2std_mat)

            sub_id_extracted = re.search('.+_subject_id_(\d+)', brain).group(1)
            if str(subject_id) in brain:
    #             print("Files for subject ",subject_id,brain,mask,atlas,tr,motion_param)
                return brain,mask,atlas,tr,motion_param,func2std_mat,MNI3mm_path

        print ('Unable to locate Subject: ',subject_id,'extracted: ',sub_id_extracted)
        # print ('Unable to locate Subject: ',subject_id)
        raise Exception('Unable to locate Subject: ',subject_id,'extracted: ',sub_id_extracted)
        # raise Exception('Unable to locate Subject: ',subject_id)
        return 0




    # Make a node
    getSubjectFilenames = Node(Function(function=get_subject_filenames, input_names=['subject_id','brain_path','mask_path','atlas_path','tr_path','motion_params_path','func2std_mat_path','MNI3mm_path'],
                                    output_names=['brain','mask','atlas','tr','motion_param','func2std_mat', 'MNI3mm_path']), name='getSubjectFilenames')


    getSubjectFilenames.inputs.brain_path = brain_path
    getSubjectFilenames.inputs.mask_path = mask_path
    getSubjectFilenames.inputs.atlas_path = atlas_path
    getSubjectFilenames.inputs.tr_path = tr_path
    getSubjectFilenames.inputs.motion_params_path = motion_params_path
    getSubjectFilenames.inputs.func2std_mat_path = func2std_mat_path
    getSubjectFilenames.inputs.MNI3mm_path = MNI3mm_path




    infosource = Node(IdentityInterface(fields=['subject_id']),
                      name="infosource")

    infosource.iterables = [('subject_id',subject_list)]



    # ## Band Pass Filtering
    # Let's do a band pass filtering on the data using the
    # code from https://neurostars.org/t/bandpass-filtering-different-outputs-from-fsl-and-nipype-custom-function/824/2

    ### AFNI

    bandpass = Node(afni.Bandpass(highpass=0.01, lowpass=0.1,
                             despike=False, no_detrend=True, notrans=True,
                             outputtype='NIFTI_GZ'),name='bandpass')

    # bandpass = Node(afni.Bandpass(highpass=0.001, lowpass=0.01,
    #                          despike=False, no_detrend=True, notrans=True,
    #                          tr=2.0,outputtype='NIFTI_GZ'),name='bandpass')


    # ## Highpass filtering

    # In[506]:

    """
    Perform temporal highpass filtering on the data
    """

    # https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dBandpass.html
    # os.chdir('/home1/varunk/Autism-Connectome-Analysis-bids-related/')

    highpass = Node(afni.Bandpass(highpass=0.009, lowpass=99999,
                             despike=False, no_detrend=True, notrans=True,
                             outputtype='NIFTI_GZ'),name='highpass')

    #  FSL bandpass/Highpass
    # highpass = Node(interface=ImageMaths(suffix='_tempfilt'),
    #                   iterfield=['in_file'],
    #                   name='highpass')
    #
    # highpass.inputs.op_string = '-bptf 27.77775001525879  -1' # 23.64 # 31.25


    # ## Smoothing
    # ### Using 6mm fwhm
    # sigma = 6/2.3548 = 2.547987090198743

    spatialSmooth = Node(interface=ImageMaths(op_string='-s 2.5479',
                                                suffix='_smoothed'),
                       name='spatialSmooth')


    # ## Performs Gram Schmidt Process
    # https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process

    # In[509]:


    def orthogonalize(in_file, mask_file):
        import numpy as np
        import nibabel as nib
        import os
        from os.path import join as opj

        def gram_schmidt(voxel_time_series, mean_vector):
            numerator = np.dot(voxel_time_series,mean_vector)
            dinominator = np.dot(mean_vector,mean_vector)
            voxel_time_series_orthogonalized = voxel_time_series - (numerator/dinominator)*mean_vector

    #         TO CONFIRM IF THE VECTORS ARE ORTHOGONAL
    #         sum_dot_prod = np.sum(np.dot(voxel_time_series_orthogonalized,mean_vector))

    #         print('Sum of entries of orthogonalized vector = ',sum_dot_prod)
            return voxel_time_series_orthogonalized


        mask_data = nib.load(mask_file)
        mask = mask_data.get_data()

        brain_data = nib.load(in_file)
        brain = brain_data.get_data()

        x_dim, y_dim, z_dim, t_dim = brain_data.shape



        # Find mean brain


        mean_vector = np.zeros(t_dim)


        num_brain_voxels = 0

        # Count the number of brain voxels
        for i in range(x_dim):
            for j in range(y_dim):
                for k in range(z_dim):
                    if mask[i,j,k] == 1:
                        mean_vector = mean_vector + brain[i,j,k,:]
                        num_brain_voxels = num_brain_voxels + 1


        mean_vector = mean_vector / num_brain_voxels

        # Orthogonalize
        for i in range(x_dim):
            for j in range(y_dim):
                for k in range(z_dim):
                    if mask[i,j,k] == 1:
                        brain[i,j,k,:] = gram_schmidt(brain[i,j,k,:], mean_vector)



        sub_id = in_file.split('/')[-1].split('.')[0].split('_')[0].split('-')[1]

        gsr_file_name = 'sub-' + sub_id + '_task-rest_run-1_bold.nii.gz'

    #     gsr_file_name_nii = gsr_file_name + '.nii.gz'

        out_file = opj(os.getcwd(),gsr_file_name) # path

        brain_with_header = nib.Nifti1Image(brain, affine=brain_data.affine,header = brain_data.header)
        nib.save(brain_with_header,gsr_file_name)

        return out_file








    # In[510]:


    globalSignalRemoval = Node(Function(function=orthogonalize, input_names=['in_file','mask_file'],
                                      output_names=['out_file']), name='globalSignalRemoval' )
    # globalSignalRemoval.inputs.mask_file = mask_file
    # globalSignalRemoval.iterables = [('in_file',file_paths)]


    # ## GLM for regression of motion parameters

    # In[511]:


    # def calc_residuals(in_file,
    #                    motion_file):
    #     """
    #     Calculates residuals of nuisance regressors -motion parameters for every voxel for a subject using GLM.
    #
    #     Parameters
    #     ----------
    #     in_file : string
    #         Path of a subject's motion corrected nifti file.
    #     motion_par_file : string
    #         path of a subject's motion parameters
    #
    #
    #     Returns
    #     -------
    #     out_file : string
    #         Path of residual file in nifti format
    #
    #     """
    #     import nibabel as nb
    #     import numpy as np
    #     import os
    #     from os.path import join as opj
    #     nii = nb.load(in_file)
    #     data = nii.get_data().astype(np.float32)
    #     global_mask = (data != 0).sum(-1) != 0
    #
    #
    #     # Check and define regressors which are provided from files
    #     if motion_file is not None:
    #         motion = np.genfromtxt(motion_file)
    #         if motion.shape[0] != data.shape[3]:
    #             raise ValueError('Motion parameters {0} do not match data '
    #                              'timepoints {1}'.format(motion.shape[0],
    #                                                      data.shape[3]))
    #         if motion.size == 0:
    #             raise ValueError('Motion signal file {0} is '
    #                              'empty'.format(motion_file))
    #
    #     # Calculate regressors
    #     regressor_map = {'constant' : np.ones((data.shape[3],1))}
    #
    #     regressor_map['motion'] = motion
    #
    #
    #     X = np.zeros((data.shape[3], 1))
    #
    #     for rname, rval in regressor_map.items():
    #         X = np.hstack((X, rval.reshape(rval.shape[0],-1)))
    #
    #     X = X[:,1:]
    #
    #     if np.isnan(X).any() or np.isnan(X).any():
    #         raise ValueError('Regressor file contains NaN')
    #
    #     Y = data[global_mask].T
    #
    #     try:
    #         B = np.linalg.inv(X.T.dot(X)).dot(X.T).dot(Y)
    #     except np.linalg.LinAlgError as e:
    #         if "Singular matrix" in e:
    #             raise Exception("Error details: {0}\n\nSingular matrix error: "
    #                             "The nuisance regression configuration you "
    #                             "selected may have been too stringent, and the "
    #                             "regression could not be completed. Ensure your "
    #                             "parameters are not too "
    #                             "extreme.\n\n".format(e))
    #         else:
    #             raise Exception("Error details: {0}\n\nSomething went wrong with "
    #                             "nuisance regression.\n\n".format(e))
    #
    #     Y_res = Y - X.dot(B)
    #
    #     data[global_mask] = Y_res.T
    #
    #     img = nb.Nifti1Image(data, header=nii.get_header(),
    #                          affine=nii.get_affine())
    #
    #     subject_name = in_file.split('/')[-1].split('.')[0]
    #     filename = subject_name + '_residual.nii.gz'
    #     out_file = os.path.join(os.getcwd(),filename )
    #     img.to_filename(out_file) # alt to nib.save
    #
    #     return out_file
    #
    #
    # # In[512]:
    #
    #
    # # Create a Node for above
    # calc_residuals = Node(Function(function=calc_residuals, input_names=['in_file','motion_file'],
    #                                 output_names=['out_file']), name='calc_residuals')

    calc_residuals = Node(Function(function=func_calc_residuals, input_names=['in_file', 'motion_file',
     'csf_mask_path', 'wm_mask_path','global_signal_flag','const','check_orthogonality'],
                                    output_names=['out_file_list']), name='calc_residuals')

    # ## Datasink
    # I needed to define the structure of what files are saved and where.

    # In[513]:


    # Create DataSink object
    dataSink = Node(DataSink(), name='datasink')

    # Name of the output folder
    dataSink.inputs.base_directory = opj(base_directory,fc_datasink_name)




    # To create the substitutions I looked the `datasink` folder where I was redirecting the output. I manually selected the part of file/folder name that I wanted to change and copied below to be substituted.
    #

    # In[514]:


    # Define substitution strings so that the data is similar to BIDS
    substitutions = [('_subject_id_', 'sub-')]

    # Feed the substitution strings to the DataSink node
    dataSink.inputs.substitutions = substitutions



    # ### Following is a Join Node that collects the preprocessed file paths and saves them in a file

    # In[516]:


    def save_file_list_function_brain(in_fc_map_brain_file):
        # Imports
        import numpy as np
        import os
        from os.path import join as opj


        file_list = np.asarray(in_fc_map_brain_file)
        print('######################## File List ######################: \n',file_list)

        np.save('fc_map_brain_file_list',file_list)
        file_name = 'fc_map_brain_file_list.npy'
        out_fc_map_brain_file = opj(os.getcwd(),file_name) # path


        return out_fc_map_brain_file



    # In[517]:


    save_file_list_brain = JoinNode(Function(function=save_file_list_function_brain, input_names=['in_fc_map_brain_file'],
                     output_names=['out_fc_map_brain_file']),
                     joinsource="infosource",
                     joinfield=['in_fc_map_brain_file'],
                     name="save_file_list_brain")

    # Utility function that saves the file paths correlation maps' (which are in npy format)

    def save_file_list_function_npy(in_fc_map_npy_file):
        # Imports
        import numpy as np
        import os
        from os.path import join as opj


        file_list = np.asarray(in_fc_map_npy_file)
        print('######################## File List ######################: \n',file_list)

        np.save('fc_map_npy_file_list',file_list)
        file_name = 'fc_map_npy_file_list.npy'
        out_fc_map_npy_file = opj(os.getcwd(),file_name) # path


        return out_fc_map_npy_file



    # In[517]:


    save_file_list_npy = JoinNode(Function(function=save_file_list_function_npy, input_names=['in_fc_map_npy_file'],
                     output_names=['out_fc_map_npy_file']),
                     joinsource="infosource",
                     joinfield=['in_fc_map_npy_file'],
                     name="save_file_list_npy")



















    # ## Create a FC node
    #
    # This node:
    # 1. Exracts the average time series of the brain ROI's using the atlas and stores
    #     it as a matrix of size [ROIs x Volumes].
    # 2. Extracts the Voxel time series and stores it in matrix of size [Voxels x Volumes]
    #


    # And save  FC matrix files in shape of brains
    def pear_coff(in_file, atlas_file, mask_file):
        # code to find how many voxels are in the brain region using the mask

            # imports
        import numpy as np
        import nibabel as nib
        import os
        from os.path import join as opj

        mask_data = nib.load(mask_file)
        mask = mask_data.get_data()

        x_dim, y_dim, z_dim = mask_data.shape


        atlasPath = atlas_file
        # Read the atlas
        atlasObject = nib.load(atlasPath)
        atlas = atlasObject.get_data()

        num_ROIs = int((np.max(atlas) - np.min(atlas) ))


        # Read the brain in_file

        brain_data = nib.load(in_file)
        brain = brain_data.get_data()

        x_dim, y_dim, z_dim, num_volumes = brain.shape


        num_brain_voxels = 0

        x_dim, y_dim, z_dim = mask_data.shape

        for i in range(x_dim):
            for j in range(y_dim):
                for k in range(z_dim):
                    if mask[i,j,k] == 1:
                        num_brain_voxels = num_brain_voxels + 1

        # Initialize a matrix of ROI time series and voxel time series

        ROI_matrix = np.zeros((num_ROIs, num_volumes))
        voxel_matrix = np.zeros((num_brain_voxels, num_volumes))

        # Fill up the voxel_matrix

        voxel_counter = 0
        for i in range(x_dim):
            for j in range(y_dim):
                for k in range(z_dim):
                    if mask[i,j,k] == 1:
                        voxel_matrix[voxel_counter,:] = brain[i,j,k,:]
                        voxel_counter = voxel_counter + 1


        # Fill up the ROI_matrix
        # Keep track of number of voxels per ROI as well by using an array - num_voxels_in_ROI[]

        num_voxels_in_ROI = np.zeros((num_ROIs,1)) # A column arrray containing number of voxels in each ROI

        for i in range(x_dim):
            for j in range(y_dim):
                for k in range(z_dim):
                    label = int(atlas[i,j,k]) - 1
                    if label != -1:
                        ROI_matrix[label,:] = np.add(ROI_matrix[label,:], brain[i,j,k,:])
                        num_voxels_in_ROI[label,0] = num_voxels_in_ROI[label,0] + 1

        ROI_matrix = np.divide(ROI_matrix,num_voxels_in_ROI) # Check if divide is working correctly

        X, Y = ROI_matrix, voxel_matrix


        # Subtract mean from X and Y

        X = np.subtract(X, np.mean(X, axis=1, keepdims=True))
        Y = np.subtract(Y, np.mean(Y, axis=1, keepdims=True))

        temp1 = np.dot(X,Y.T)
        temp2 = np.sqrt(np.sum(np.multiply(X,X), axis=1, keepdims=True))
        temp3 = np.sqrt(np.sum(np.multiply(Y,Y), axis=1, keepdims=True))
        temp4 = np.dot(temp2,temp3.T)
        coff_matrix = np.divide(temp1, (temp4 + 1e-7))


        # Check if any ROI is missing and replace the NAN values in coff_matrix by 0
        if np.argwhere(np.isnan(coff_matrix)).shape[0] != 0:
            print("Some ROIs are not present. Replacing NAN in coff matrix by 0")
            np.nan_to_num(coff_matrix, copy=False)

        # TODO: when I have added 1e-7 in the dinominator, then why did I feel the need to replace NAN by zeros
        sub_id = in_file.split('/')[-1].split('.')[0].split('_')[0].split('-')[1]


        fc_file_name = sub_id + '_fc_map'

        print ("Pear Matrix calculated for subject: ",sub_id)

        roi_brain_matrix = coff_matrix
        brain_file = in_file


        x_dim, y_dim, z_dim, t_dim = brain.shape

        (brain_data.header).set_data_shape([x_dim,y_dim,z_dim,num_ROIs])

        brain_roi_tensor = np.zeros((brain_data.header.get_data_shape()))

        print("Creating brain for Subject-",sub_id)
        for roi in range(num_ROIs):
            brain_voxel_counter = 0
            for i in range(x_dim):
                for j in range(y_dim):
                    for k in range(z_dim):
                        if mask[i,j,k] == 1:
                            brain_roi_tensor[i,j,k,roi] = roi_brain_matrix[roi,brain_voxel_counter]
                            brain_voxel_counter = brain_voxel_counter + 1


            assert (brain_voxel_counter == len(roi_brain_matrix[roi,:]))
        print("Created brain for Subject-",sub_id)


        path = os.getcwd()
        fc_brain_file_name = fc_file_name + '.nii.gz'
        out_file_brain = opj(path,fc_brain_file_name)

        brain_with_header = nib.Nifti1Image(brain_roi_tensor, affine=brain_data.affine,header = brain_data.header)
        nib.save(brain_with_header,out_file_brain)


        fc_map_brain_file = out_file_brain

        return fc_map_brain_file



    # In[521]:


    # Again Create the Node and set default values to paths

    pearcoff = Node(Function(function=pear_coff, input_names=['in_file','atlas_file','mask_file'],
                                    output_names=['fc_map_brain_file']), name='pearcoff')



    # # IMPORTANT:
    # * The ROI 255 has been removed due to resampling. Therefore the FC maps will have nan at that row. So don't use that ROI :)
    # * I came to know coz I keep getting this error: RuntimeWarning: invalid value encountered in true_divide
    # * To debug it, I read the coff matrix and checked its diagnol to discover the nan value.
    #
    #
    #

    # ## Extract volumes




    # ExtractROI - For volCorrect
    extract = Node(ExtractROI(t_size=-1),
                   output_type='NIFTI',
                   name="extract")



    # ###  Node for applying xformation matrix to functional data
    #

    # In[523]:


    func2std_xform = Node(FLIRT(output_type='NIFTI_GZ',
                             apply_xfm=True), name="func2std_xform")



    def convert_nii_2_npy_func(in_file):
        import nibabel as  nib
        import numpy as np
        import os
        from os.path import join as opj

        file_name_nii = in_file
        subject_id = file_name_nii.split('/')[-1].split('_')[0]
        filename = 'sub-'+subject_id+'.npy'
        print('converting the Brain of subject %s to npy with filename %s'%(subject_id, filename))

        brain = nib.load(file_name_nii).get_data()
        np.save(filename,brain)


        path = os.getcwd()
        out_file = opj(path,filename)

        return out_file


    nii2npy = Node(Function(function=convert_nii_2_npy_func, input_names=['in_file'],
                                    output_names=['out_file']), name='nii2npy')







    # motion_param_regression = 1
    # band_pass_filtering = 0
    # global_signal_regression = 0
    # smoothing = 1
    # volcorrect = 1
    if num_proc == None:
        num_proc = 7

    combination = 'motionRegress' + str(int(motion_param_regression)) + \
     'global' + str(int(global_signal_regression)) + 'smoothing' + str(int(smoothing)) +\
     'filt' + str(int(band_pass_filtering))

    print("Combination: ",combination)

    binary_string = str(int(motion_param_regression)) + str(int(global_signal_regression)) + \
    str(int(smoothing)) + str(int(band_pass_filtering)) + str(int(volcorrect))

    base_dir = opj(base_directory,functional_connectivity_directory)
    # wf = Workflow(name=functional_connectivity_directory)
    wf = Workflow(name=combination)

    wf.base_dir = base_dir # Dir where all the outputs will be stored.

    wf.connect(infosource ,'subject_id', getSubjectFilenames, 'subject_id')


    # ------- Dynamic Pipeline ------------------------


    nodes = [
    calc_residuals,
    globalSignalRemoval,
    spatialSmooth,
    bandpass,
    volCorrect]


    # from nipype.interfaces import fsl

    old_node = getSubjectFilenames
    old_node_output = 'brain'

    binary_string = binary_string+'0' # so that the loop runs one more time
    for idx, include in enumerate(binary_string):
        # 11111
        # motion_param_regression
        # global_signal_regression
        # smoothing
        # band_pass_filtering
        # volcorrect

        if old_node == calc_residuals:
            old_node_output = 'out_file_list'
        elif old_node == extract :
            old_node_output = 'roi_file'
        elif old_node == globalSignalRemoval:
            old_node_output = 'out_file'
        elif old_node == bandpass:
            old_node_output = 'out_file'
        elif old_node == highpass:
            old_node_output = 'out_file'
        elif old_node == spatialSmooth:
            old_node_output = 'out_file'
        elif old_node == volCorrect:
            old_node_output = 'out_file'


        if int(include):
            # if old_node is None:
            #
            #     wf.add_nodes([nodes[idx]])
            #
            # else:



            new_node = nodes[idx]


            if new_node == calc_residuals:
                wf.connect([(getSubjectFilenames, calc_residuals, [('motion_param', 'motion_file')])])
                new_node_input = 'in_file'

            elif new_node == extract :
                wf.connect([( volCorrect, extract, [('t_min','t_min')])])
                new_node_input = 'in_file'

            elif new_node == globalSignalRemoval:
                wf.connect([(getSubjectFilenames, globalSignalRemoval, [('mask','mask_file')])])
                new_node_input = 'in_file'

            elif new_node == bandpass:
                wf.connect([(getSubjectFilenames, bandpass, [('tr','tr')])])
                new_node_input = 'in_file'

            elif new_node == highpass:
                wf.connect([(getSubjectFilenames, highpass, [('tr','tr')])]) #Commenting for FSL
                new_node_input = 'in_file'

            elif new_node == spatialSmooth:
                new_node_input = 'in_file'

            elif new_node == volCorrect:
                wf.connect([(infosource, volCorrect, [('subject_id','sub_id')])])
                wf.connect([( volCorrect, extract, [('t_min','t_min')])])
                new_node = extract
                new_node_input = 'in_file'


            wf.connect(old_node, old_node_output, new_node, new_node_input)

            old_node = new_node


        else:
            if idx == 3: # bandpas == 0 => Highpass
                new_node = highpass
                wf.connect([(getSubjectFilenames, highpass, [('tr','tr')])]) #Commenting for FSL
                new_node_input = 'in_file'

                wf.connect(old_node, old_node_output, new_node, new_node_input)

                old_node = new_node

    wf.connect(old_node, old_node_output, pearcoff, 'in_file')
    wf.connect(getSubjectFilenames,'atlas', pearcoff, 'atlas_file')
    wf.connect(getSubjectFilenames, 'mask', pearcoff, 'mask_file')

    wf.connect(pearcoff, 'fc_map_brain_file', func2std_xform ,'in_file')
    wf.connect(getSubjectFilenames,'func2std_mat', func2std_xform, 'in_matrix_file')
    wf.connect(getSubjectFilenames, 'MNI3mm_path', func2std_xform,'reference')

    folder_name = combination + '.@fc_map_brain_file'
    wf.connect(func2std_xform, 'out_file',  save_file_list_brain, 'in_fc_map_brain_file')
    wf.connect(save_file_list_brain, 'out_fc_map_brain_file',  dataSink,folder_name)

    if save_npy == 1:
        folder_name = combination + '.@fc_map_npy_file'
        wf.connect(func2std_xform, 'out_file', nii2npy,'in_file')
        wf.connect(nii2npy, 'out_file',  save_file_list_npy, 'in_fc_map_npy_file')
        wf.connect(save_file_list_npy, 'out_fc_map_npy_file',  dataSink,folder_name)



    TEMP_DIR_FOR_STORAGE = opj(base_directory,'crash_files')
    wf.config = {"execution": {"crashdump_dir": TEMP_DIR_FOR_STORAGE}}

    wf.write_graph(graph2use='flat', format='png')
    wf.run('MultiProc', plugin_args={'n_procs': num_proc})

    # -------------------------------------------------
