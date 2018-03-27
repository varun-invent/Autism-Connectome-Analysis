
# coding: utf-8

# ## Post processing *Version 3*:
# * Global Signal Regression using orthogonalization
# * Band Pass filtering 0.1 - 0.01 Hz
# * Motion regression using GLM
#
#

# In[752]:


from bids.grabbids import BIDSLayout
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


# In[753]:


# Paths
def main(motion_param_regression=0, global_signal_regression=0, band_pass_filtering=0, smoothing=0,  number_of_subjects = -1, functional_connectivity_directory =None):
    path_cwd = os.getcwd()
    path_split_list = path_cwd.split('/')
    s = path_split_list[0:-1] # for getting to the parent dir of pwd
    s = opj('/',*s) # *s converts list to path, # very important to add '/' in the begining so it is read as directory later


    num_proc = 7

    # In[754]:



    # json_path = opj(data_directory,'task-rest_bold.json')

    json_path = 'scripts/json/paths.json'
    with open(json_path, 'rt') as fp:
        task_info = json.load(fp)



    # In[755]:


    # base_directory = opj(s,'result')
    # parent_wf_directory = 'preprocessPipeline_ABIDE2_GU1_withfloat'
    # child_wf_directory = 'coregistrationPipeline'

    # data_directory = opj(s,"data/ABIDE2-BIDS/GU1")

    # datasink_name = 'datasink_preprocessed_ABIDE2_GU1_withfloat'
    if functional_connectivity_directory == None:
        functional_connectivity_directory = task_info["functional_connectivity_directory"]


    base_directory = opj(s,task_info["base_directory_for_results"])
    motion_correction_bet_directory = task_info["motion_correction_bet_directory"]
    parent_wf_directory = task_info["parent_wf_directory"]
    coreg_reg_directory = task_info["coreg_reg_directory"]
    atlas_resize_reg_directory = task_info["atlas_resize_reg_directory"]
    data_directory = opj(s,task_info["data_directory"])
    datasink_name = task_info["datasink_name"]
    fc_datasink_name = task_info["fc_datasink_name"]

    atlasPath = opj(s,task_info["atlas_path"])

    # mask_file = '/media/varun/LENOVO4/Projects/result/preprocessPipeline/coregistrationPipeline/_subject_id_0050952/skullStrip/sub-0050952_T1w_resample_brain_mask.nii.gz'
    # os.chdir(path)


    # In[756]:


    brain_path = opj(base_directory,datasink_name,'preprocessed_brain_paths/brain_file_list.npy')
    mask_path = opj(base_directory,datasink_name,'preprocessed_mask_paths/mask_file_list.npy')
    atlas_path = opj(base_directory,datasink_name,'atlas_paths/atlas_file_list.npy')
    tr_path = opj(base_directory,datasink_name,'tr_paths/tr_list.npy')
    motion_params_path = opj(base_directory,datasink_name,'motion_params_paths/motion_params_file_list.npy')
    func2std_mat_path = opj(base_directory, datasink_name,'joint_xformation_matrix_paths/joint_xformation_matrix_file_list.npy')
    MNI3mm_path = opj(base_directory,parent_wf_directory,motion_correction_bet_directory,coreg_reg_directory,'resample_mni/MNI152_T1_2mm_brain_resample.nii')



    # brain_list = np.load('../results_again_again/ABIDE1_Preprocess_Datasink/preprocessed_brain_paths/brain_file_list.npy')


    # In[757]:


    # brain_path,mask_path,atlas_path,tr_path,motion_params_path


    # In[758]:


    brain_path = np.load(brain_path)
    mask_path = np.load(mask_path)
    atlas_path = np.load(atlas_path)
    tr_path = np.load(tr_path)
    motion_params_path = np.load(motion_params_path)

    func2std_mat_path = np.load(func2std_mat_path)

    # In[759]:


    # for a,b,c,d,e in zip(brain_path,mask_path,atlas_path,tr_path,motion_params_path):
    #     print (a,b,c,d,e,'\n')


    # In[760]:


    layout = BIDSLayout(data_directory)

    # number_of_subjects = 1 # Number of subjects you wish to preprocess

    if number_of_subjects ==-1:
        number_of_subjects = len(layout.get_subjects())

    print('Number of Subjects: ',number_of_subjects)


    # Checking the Data directory Structure

    # In[761]:


    # len(layout.get_subjects()) # working!Gives us list of all the subjects


    # In[762]:


    # layout.get_subjects();


    # To get the metadata associated with a subject. [Takes as argumment the filename of subject ]

    # Create a list of subjects

    # In[763]:


    subject_list = (layout.get_subjects())[0:number_of_subjects]


    # In[764]:


    # layout.get();


    # Create our own custom function - BIDSDataGrabber using a Function Interface.

    # In[765]:


    # def get_nifti_filenames(subject_id,data_dir):
    # #     Remember that all the necesary imports need to be INSIDE the function for the Function Interface to work!
    #     from bids.grabbids import BIDSLayout
    #
    #     layout = BIDSLayout(data_dir)
    #
    #     anat_file_path = [f.filename for f in layout.get(subject=subject_id, type='T1w', extensions=['nii', 'nii.gz'])]
    #     func_file_path = [f.filename for f in layout.get(subject=subject_id, type='bold', run='1', extensions=['nii', 'nii.gz'])]
    #
    #     return anat_file_path[0],func_file_path[0]
    #
    # # Refer to Supplementary material section One for info on arguments for layout.get()
    #
    #
    # # Wrap it inside a Node
    #
    # # In[766]:
    #
    #
    # BIDSDataGrabber = Node(Function(function=get_nifti_filenames, input_names=['subject_id','data_dir'],
    #                                 output_names=['anat_file_path','func_file_path']), name='BIDSDataGrabber')
    # # BIDSDataGrabber.iterables = [('subject_id',subject_list)]
    # BIDSDataGrabber.inputs.data_dir = data_directory
    #
    #
    # # In[767]:
    #
    #
    # len(subject_list)


    # ## Define a function to fetch the filenames of a particular subject ID

    # In[768]:


    def get_subject_filenames(subject_id,brain_path,mask_path,atlas_path,tr_path,motion_params_path,func2std_mat_path,MNI3mm_path):
        import re

        for brain,mask,atlas,tr,motion_param,func2std_mat in zip(brain_path,mask_path,atlas_path,tr_path,motion_params_path,func2std_mat_path):
            sub_id_extracted = re.search('.+_subject_id_(\d+)', brain).group(1)
            if subject_id == sub_id_extracted:
    #             print("Files for subject ",subject_id,brain,mask,atlas,tr,motion_param)
                return brain,mask,atlas,tr,motion_param,func2std_mat,MNI3mm_path

        print ('Unable to locate Subject: ',subject_id,'extracted: ',sub_id_extracted)
        return 0



    # In[769]:


    # Make a node
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

    # In[770]:


    # import re
    # text = '/home1/varunk/results_again_again/ABIDE1_Preprocess/motion_correction_bet/coreg_reg/atlas_resize_reg_directory/_subject_id_0050004/111std2func_xform/fullbrain_atlas_thr0-2mm_resample_flirt.nii'

    # try:
    #     found = re.search('.+_subject_id_(\d+)', text).group(1)
    # except AttributeError:
    #     # AAA, ZZZ not found in the original string
    #     found = '' # apply your error handling

    # # found: 1234


    # In[771]:


    # found


    # In[772]:


    infosource = Node(IdentityInterface(fields=['subject_id']),
                      name="infosource")

    infosource.iterables = [('subject_id',subject_list)]


    # ,'brain_path','mask_path','atlas_path','tr_path','motion_params_path'
    # infosource.brain_path = brain_path
    # infosource.mask_path = mask_path
    # infosource.atlas_path = atlas_path
    # infosource.tr_path = tr_path
    # infosource.motion_params_path = motion_params_path


    # ## Band Pass Filtering
    # Let's do a band pass filtering on the data using the code from https://neurostars.org/t/bandpass-filtering-different-outputs-from-fsl-and-nipype-custom-function/824/2

    # In[773]:


    ### AFNI

    bandpass = Node(afni.Bandpass(highpass=0.01, lowpass=0.1,
                             despike=False, no_detrend=True, notrans=True,
                             outputtype='NIFTI_GZ'),name='bandpass')

    # bandpass = Node(afni.Bandpass(highpass=0.001, lowpass=0.01,
    #                          despike=False, no_detrend=True, notrans=True,
    #                          tr=2.0,outputtype='NIFTI_GZ'),name='bandpass')


    # bandpass.inputs.mask = MNI152_2mm.outputs.mask_file

    # highpass=0.008, lowpass=0.08,


    # In[774]:


    # Testing bandpass on the func data in subject's space

    # First comment out the bandpass.inputs.mask as it is in standard space.

    # subject_id = layout.get_subjects()[0] # gives the first subject's ID
    # func_file_path = [f.filename for f in layout.get(subject=subject_id, type='bold', extensions=['nii', 'nii.gz'])]
    # bandpass.inputs.in_file = func_file_path[0]
    # res = bandpass.run();


    # In[775]:
    # Highpass filter
    highpass = Node(afni.Bandpass(highpass=0.01, lowpass=99999,
                         despike=False, no_detrend=True, notrans=True,
                         outputtype='NIFTI_GZ'),name='highpass')

    # # Test

    # subject_id = layout.get_subjects()[0] # gives the first subject's ID
    # func_file_path = [f.filename for f in layout.get(subject=subject_id, type='bold', extensions=['nii', 'nii.gz'])]
    # highpass.inputs.in_file = func_file_path[0]
    # res = highpass.run();

    # res.outputs.out_file


    # ## Performs Gram Schmidt Process
    # https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process

    # In[776]:


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








    # In[777]:


    globalSignalRemoval = Node(Function(function=orthogonalize, input_names=['in_file','mask_file'],
                                      output_names=['out_file']), name='globalSignalRemoval' )
    # globalSignalRemoval.inputs.mask_file = mask_file
    # globalSignalRemoval.iterables = [('in_file',file_paths)]


    # ## GLM for regression of motion parameters

    # In[778]:


    def calc_residuals(subject,
                       motion_file):
        """
        Calculates residuals of nuisance regressors -motion parameters for every voxel for a subject using GLM.

        Parameters
        ----------
        subject : string
            Path of a subject's motion corrected nifti file.
        motion_par_file : string
            path of a subject's motion parameters


        Returns
        -------
        residual_file : string
            Path of residual file in nifti format

        """
        import nibabel as nb
        import numpy as np
        import os
        from os.path import join as opj
        nii = nb.load(subject)
        data = nii.get_data().astype(np.float32)
        global_mask = (data != 0).sum(-1) != 0


        # Check and define regressors which are provided from files
        if motion_file is not None:
            motion = np.genfromtxt(motion_file)
            if motion.shape[0] != data.shape[3]:
                raise ValueError('Motion parameters {0} do not match data '
                                 'timepoints {1}'.format(motion.shape[0],
                                                         data.shape[3]))
            if motion.size == 0:
                raise ValueError('Motion signal file {0} is '
                                 'empty'.format(motion_file))

        # Calculate regressors
        regressor_map = {'constant' : np.ones((data.shape[3],1))}

        regressor_map['motion'] = motion


        X = np.zeros((data.shape[3], 1))

        for rname, rval in regressor_map.items():
            X = np.hstack((X, rval.reshape(rval.shape[0],-1)))

        X = X[:,1:]

        if np.isnan(X).any() or np.isnan(X).any():
            raise ValueError('Regressor file contains NaN')

        Y = data[global_mask].T

        try:
            B = np.linalg.inv(X.T.dot(X)).dot(X.T).dot(Y)
        except np.linalg.LinAlgError as e:
            if "Singular matrix" in e:
                raise Exception("Error details: {0}\n\nSingular matrix error: "
                                "The nuisance regression configuration you "
                                "selected may have been too stringent, and the "
                                "regression could not be completed. Ensure your "
                                "parameters are not too "
                                "extreme.\n\n".format(e))
            else:
                raise Exception("Error details: {0}\n\nSomething went wrong with "
                                "nuisance regression.\n\n".format(e))

        Y_res = Y - X.dot(B)

        data[global_mask] = Y_res.T

        img = nb.Nifti1Image(data, header=nii.get_header(),
                             affine=nii.get_affine())

        subject_name = subject.split('/')[-1].split('.')[0]
        filename = subject_name + '_residual.nii.gz'
        residual_file = os.path.join(os.getcwd(),filename )
        img.to_filename(residual_file) # alt to nib.save

        return residual_file


    # In[779]:


    # Create a Node for above
    calc_residuals = Node(Function(function=calc_residuals, input_names=['subject','motion_file'],
                                    output_names=['residual_file']), name='calc_residuals')


    # ## Datasink
    # I needed to define the structure of what files are saved and where.

    # In[780]:


    # Create DataSink object
    dataSink = Node(DataSink(), name='datasink')

    # Name of the output folder
    dataSink.inputs.base_directory = opj(base_directory,fc_datasink_name)




    # To create the substitutions I looked the `datasink` folder where I was redirecting the output. I manually selected the part of file/folder name that I wanted to change and copied below to be substituted.
    #

    # In[781]:


    # Define substitution strings so that the data is similar to BIDS
    substitutions = [('_subject_id_', 'sub-')]

    # Feed the substitution strings to the DataSink node
    dataSink.inputs.substitutions = substitutions

    # ('_resample_brain_flirt.nii_brain', ''),
    # ('_roi_st_mcf_flirt.nii_brain_flirt', ''),



    # In[782]:


    # base_directory


    # ### Following is a Join Node that collects the preprocessed file paths and saves them in a file

    # In[783]:


    def save_file_list_function(in_fc_map_brain_file):
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



    # In[784]:


    save_file_list = JoinNode(Function(function=save_file_list_function, input_names=['in_fc_map_brain_file'],
                     output_names=['out_fc_map_brain_file']),
                     joinsource="infosource",
                     joinfield=['in_fc_map_brain_file'],
                     name="save_file_list")


    # ## Create a FC node
    #
    # This node:
    # 1. Exracts the average time series of the brain ROI's using the atlas and stores
    #     it as a matrix of size [ROIs x Volumes].
    # 2. Extracts the Voxel time series and stores it in matrix of size [Voxels x Volumes]
    #

    # In[785]:


    # Saves the brains instead of FC matrix files
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
        fc_file_name = fc_file_name + '.nii.gz'
        out_file = opj(path,fc_file_name)

        brain_with_header = nib.Nifti1Image(brain_roi_tensor, affine=brain_data.affine,header = brain_data.header)
        nib.save(brain_with_header,out_file)


        fc_map_brain_file = out_file
        return fc_map_brain_file





    # In[786]:


    # Again Create the Node and set default values to paths

    pearcoff = Node(Function(function=pear_coff, input_names=['in_file','atlas_file','mask_file'],
                                    output_names=['fc_map_brain_file']), name='pearcoff')

    # pearcoff.inputs.atlas_file = atlasPath
    # pearcoff.inputs.num_brain_voxels = num_brain_voxels
    # pearcoff.inputs.mask_file = mask_file


    # # IMPORTANT:
    # * The ROI 255 has been removed due to resampling. Therefore the FC maps will have nan at that row. So don't use that ROI :)
    # * I came to know coz I keep getting this error: RuntimeWarning: invalid value encountered in true_divide
    # * To debug it, I read the coff matrix and checked its diagnol to discover the nan value.
    #
    #
    #

    # In[787]:


    # %%time
    # pearcoff.run()
    # Spatial smoothing of 6 fwhm .i.e. sigma = 2.5479
    spatialSmooth = Node(interface=ImageMaths(op_string='-s 2.5479',
                                            suffix='_smoothed'),
                   name='spatialSmooth')

    # In[794]:

    #
    # motion_param_reg = [True, False]
    # global_signal_reg = [True, False]
    # band_pass_filt= [True, False]
    # for motion_param_regression, global_signal_regression, band_pass_filtering in zip(motion_param_reg, global_signal_reg, band_pass_filt):
    #     print (motion_param_regression, global_signal_regression, band_pass_filtering)


    # In[801]:

    #
    # import itertools
    # itr = (list(itertools.product([0, 1], repeat=3)))
    #
    # for motion_param_regression, global_signal_regression, band_pass_filtering in itr:
    #     print(motion_param_regression, global_signal_regression, band_pass_filtering)


    # In[802]:


    # import itertools
    # itr = (list(itertools.product([0, 1], repeat=3)))
    #
    # for motion_param_regression, global_signal_regression, band_pass_filtering in itr:

    # For transforming FUnctional to Standard space
    func2std_xform = Node(FLIRT(output_type='NIFTI_GZ',
                         apply_xfm=True), name="func2std_xform")


    # smoothing = 1

    combination = 'motionRegress' + str(int(motion_param_regression)) + 'filt' + \
              str(int(band_pass_filtering)) + 'global' + str(int(global_signal_regression)) + \
              'smoothing' + str(int(smoothing))

    print("Combination: ",combination)

    # new_base_dir = opj(base_directory, functional_connectivity_directory)
    wf = Workflow(name=functional_connectivity_directory)
    # base_dir = opj(s,'result')
    wf.base_dir = base_directory # Dir where all the outputs will be stored(inside BETFlow folder).

    wf.connect([(infosource , getSubjectFilenames, [('subject_id','subject_id')])])

    if motion_param_regression == 1 and global_signal_regression == 1 and band_pass_filtering == 1 and smoothing == 1: #111
        wf.connect([(getSubjectFilenames, calc_residuals, [('brain','subject')])])
        wf.connect([(getSubjectFilenames, calc_residuals, [('motion_param', 'motion_file')])])

        wf.connect([(calc_residuals, globalSignalRemoval, [('residual_file','in_file')] )])
        wf.connect([(getSubjectFilenames, globalSignalRemoval, [('mask','mask_file')])])

        wf.connect([(globalSignalRemoval, bandpass, [('out_file','in_file')])])
        wf.connect([(getSubjectFilenames, bandpass, [('tr','tr')])])

        wf.connect([( bandpass, spatialSmooth, [('out_file','in_file')])])

        wf.connect([( spatialSmooth, pearcoff, [('out_file','in_file')])])

        # wf.connect([( bandpass, pearcoff, [('out_file','in_file')])])
        wf.connect([( getSubjectFilenames, pearcoff, [('atlas','atlas_file')])])
        wf.connect([( getSubjectFilenames, pearcoff, [('mask','mask_file')])])

        # ---------------------------------------------------------------------------------------
        wf.connect([(pearcoff, func2std_xform, [('fc_map_brain_file','in_file')])])
        wf.connect([(getSubjectFilenames, func2std_xform, [('func2std_mat','in_matrix_file')])])
        wf.connect([(getSubjectFilenames, func2std_xform, [('MNI3mm_path','reference')])])

#         -- send out file to save file list and then save the outputs



        folder_name = 'pearcoff_' + combination + '.@fc_map_brain_file'



        wf.connect([(func2std_xform,  save_file_list, [('out_file','in_fc_map_brain_file')])])
        # --------------------------------------------------------------------------------------------


        wf.connect([(save_file_list,  dataSink, [('out_fc_map_brain_file',folder_name)])])


        #  wf.connect([(bandpass,  dataSink, [('out_file','motionRegress_filt_global.@out_file')])])


        # if motion_param_regression == 1 and global_signal_regression == 1:
        wf.write_graph(graph2use='flat', format='png')
        wf.run('MultiProc', plugin_args={'n_procs': num_proc})


    elif motion_param_regression == 1 and global_signal_regression == 1 and band_pass_filtering == 0 and smoothing == 1: # 110
            wf.connect([(getSubjectFilenames, calc_residuals, [('brain','subject')])])
            wf.connect([(getSubjectFilenames, calc_residuals, [('motion_param', 'motion_file')])])

            wf.connect([(calc_residuals, globalSignalRemoval, [('residual_file','in_file')] )])
            wf.connect([(getSubjectFilenames, globalSignalRemoval, [('mask','mask_file')])])

            wf.connect([(globalSignalRemoval, highpass, [('out_file','in_file')])])
            wf.connect([(getSubjectFilenames, highpass, [('tr','tr')])])

            wf.connect([( bandpass, spatialSmooth, [('out_file','in_file')])])

            wf.connect([( spatialSmooth, pearcoff, [('out_file','in_file')])])


            # wf.connect([( highpass, pearcoff, [('out_file','in_file')])])
            wf.connect([( getSubjectFilenames, pearcoff, [('atlas','atlas_file')])])
            wf.connect([( getSubjectFilenames, pearcoff, [('mask','mask_file')])])

            # ---------------------------------------------------------------------------------------
            wf.connect([(pearcoff, func2std_xform, [('fc_map_brain_file','in_file')])])
            wf.connect([(getSubjectFilenames, func2std_xform, [('func2std_mat','in_matrix_file')])])
            wf.connect([(getSubjectFilenames, func2std_xform, [('MNI3mm_path','reference')])])

    #         -- send out file to save file list and then save the outputs



            folder_name = 'pearcoff_' + combination + '.@fc_map_brain_file'



            wf.connect([(func2std_xform,  save_file_list, [('out_file','in_fc_map_brain_file')])])
            # --------------------------------------------------------------------------------------------


            wf.connect([(save_file_list,  dataSink, [('out_fc_map_brain_file',folder_name)])])

            wf.write_graph(graph2use='flat', format='png')
            wf.run('MultiProc', plugin_args={'n_procs': num_proc})


    elif motion_param_regression == 1 and global_signal_regression == 0 and band_pass_filtering == 1 and smoothing == 1: # 101
            wf.connect([(getSubjectFilenames, calc_residuals, [('brain','subject')])])
            wf.connect([(getSubjectFilenames, calc_residuals, [('motion_param', 'motion_file')])])

    #         wf.connect([(calc_residuals, globalSignalRemoval, [('residual_file','in_file')] )])
    #         wf.connect([(getSubjectFilenames, globalSignalRemoval, [('mask','mask_file')])])

            wf.connect([(calc_residuals, bandpass, [('residual_file','in_file')])])
            wf.connect([(getSubjectFilenames, bandpass, [('tr','tr')])])

            wf.connect([( bandpass, spatialSmooth, [('out_file','in_file')])])

            wf.connect([( spatialSmooth, pearcoff, [('out_file','in_file')])])


            # wf.connect([( bandpass, pearcoff, [('out_file','in_file')])])
            wf.connect([( getSubjectFilenames, pearcoff, [('atlas','atlas_file')])])
            wf.connect([( getSubjectFilenames, pearcoff, [('mask','mask_file')])])

            # ---------------------------------------------------------------------------------------
            wf.connect([(pearcoff, func2std_xform, [('fc_map_brain_file','in_file')])])
            wf.connect([(getSubjectFilenames, func2std_xform, [('func2std_mat','in_matrix_file')])])
            wf.connect([(getSubjectFilenames, func2std_xform, [('MNI3mm_path','reference')])])

    #         -- send out file to save file list and then save the outputs



            folder_name = 'pearcoff_' + combination + '.@fc_map_brain_file'



            wf.connect([(func2std_xform,  save_file_list, [('out_file','in_fc_map_brain_file')])])
            # --------------------------------------------------------------------------------------------


            wf.connect([(save_file_list,  dataSink, [('out_fc_map_brain_file',folder_name)])])

            wf.write_graph(graph2use='flat', format='png')
            wf.run('MultiProc', plugin_args={'n_procs': num_proc})

    elif motion_param_regression == 1 and global_signal_regression == 0 and band_pass_filtering == 0 and smoothing == 1: # 100
            wf.connect([(getSubjectFilenames, calc_residuals, [('brain','subject')])])
            wf.connect([(getSubjectFilenames, calc_residuals, [('motion_param', 'motion_file')])])

    #         wf.connect([(calc_residuals, globalSignalRemoval, [('residual_file','in_file')] )])
    #         wf.connect([(getSubjectFilenames, globalSignalRemoval, [('mask','mask_file')])])

            wf.connect([(calc_residuals, highpass, [('residual_file','in_file')])])
            wf.connect([(getSubjectFilenames, highpass, [('tr','tr')])])

            wf.connect([( highpass, spatialSmooth, [('out_file','in_file')])])
            wf.connect([( spatialSmooth, pearcoff, [('out_file','in_file')])])


            # wf.connect([( highpass, pearcoff, [('out_file','in_file')])])
            wf.connect([( getSubjectFilenames, pearcoff, [('atlas','atlas_file')])])
            wf.connect([( getSubjectFilenames, pearcoff, [('mask','mask_file')])])

            # ---------------------------------------------------------------------------------------
            wf.connect([(pearcoff, func2std_xform, [('fc_map_brain_file','in_file')])])
            wf.connect([(getSubjectFilenames, func2std_xform, [('func2std_mat','in_matrix_file')])])
            wf.connect([(getSubjectFilenames, func2std_xform, [('MNI3mm_path','reference')])])

    #         -- send out file to save file list and then save the outputs



            folder_name = 'pearcoff_' + combination + '.@fc_map_brain_file'



            wf.connect([(func2std_xform,  save_file_list, [('out_file','in_fc_map_brain_file')])])
            # --------------------------------------------------------------------------------------------

            wf.connect([(save_file_list,  dataSink, [('out_fc_map_brain_file',folder_name)])])

            wf.write_graph(graph2use='flat', format='png')
            wf.run('MultiProc', plugin_args={'n_procs': num_proc})


    elif motion_param_regression == 0 and global_signal_regression == 1 and band_pass_filtering == 1: # 011
    #         wf.connect([(getSubjectFilenames, calc_residuals, [('brain','subject')])])
    #         wf.connect([(getSubjectFilenames, calc_residuals, [('motion_param', 'motion_file')])])

            wf.connect([(getSubjectFilenames, globalSignalRemoval, [('brain','in_file')] )])
            wf.connect([(getSubjectFilenames, globalSignalRemoval, [('mask','mask_file')])])

            wf.connect([(globalSignalRemoval, bandpass, [('out_file','in_file')])])
            wf.connect([(getSubjectFilenames, bandpass, [('tr','tr')])])

            wf.connect([( bandpass, pearcoff, [('out_file','in_file')])])
            wf.connect([( getSubjectFilenames, pearcoff, [('atlas','atlas_file')])])
            wf.connect([( getSubjectFilenames, pearcoff, [('mask','mask_file')])])


            # ---------------------------------------------------------------------------------------
            wf.connect([(pearcoff, func2std_xform, [('fc_map_brain_file','in_file')])])
            wf.connect([(getSubjectFilenames, func2std_xform, [('func2std_mat','in_matrix_file')])])
            wf.connect([(getSubjectFilenames, func2std_xform, [('MNI3mm_path','reference')])])

    #         -- send out file to save file list and then save the outputs



            folder_name = 'pearcoff_' + combination + '.@fc_map_brain_file'



            wf.connect([(func2std_xform,  save_file_list, [('out_file','in_fc_map_brain_file')])])
            # --------------------------------------------------------------------------------------------

            wf.connect([(save_file_list,  dataSink, [('out_fc_map_brain_file',folder_name)])])

            wf.write_graph(graph2use='flat', format='png')
            wf.run('MultiProc', plugin_args={'n_procs': num_proc})



    elif motion_param_regression == 0 and global_signal_regression == 1 and band_pass_filtering == 0: # 010
    #         wf.connect([(getSubjectFilenames, calc_residuals, [('brain','subject')])])
    #         wf.connect([(getSubjectFilenames, calc_residuals, [('motion_param', 'motion_file')])])

            wf.connect([(getSubjectFilenames, globalSignalRemoval, [('brain','in_file')] )])
            wf.connect([(getSubjectFilenames, globalSignalRemoval, [('mask','mask_file')])])

            wf.connect([(globalSignalRemoval, highpass, [('out_file','in_file')])])
            wf.connect([(getSubjectFilenames, highpass, [('tr','tr')])])

            wf.connect([( highpass, pearcoff, [('out_file','in_file')])])
            wf.connect([( getSubjectFilenames, pearcoff, [('atlas','atlas_file')])])
            wf.connect([( getSubjectFilenames, pearcoff, [('mask','mask_file')])])


            # ---------------------------------------------------------------------------------------
            wf.connect([(pearcoff, func2std_xform, [('fc_map_brain_file','in_file')])])
            wf.connect([(getSubjectFilenames, func2std_xform, [('func2std_mat','in_matrix_file')])])
            wf.connect([(getSubjectFilenames, func2std_xform, [('MNI3mm_path','reference')])])

    #         -- send out file to save file list and then save the outputs



            folder_name = 'pearcoff_' + combination + '.@fc_map_brain_file'



            wf.connect([(func2std_xform,  save_file_list, [('out_file','in_fc_map_brain_file')])])
            # --------------------------------------------------------------------------------------------

            wf.connect([(save_file_list,  dataSink, [('out_fc_map_brain_file',folder_name)])])

            wf.write_graph(graph2use='flat', format='png')
            wf.run('MultiProc', plugin_args={'n_procs': num_proc})


    elif motion_param_regression == 0 and global_signal_regression == 0 and band_pass_filtering == 1: # 001
    #         wf.connect([(getSubjectFilenames, calc_residuals, [('brain','subject')])])
    #         wf.connect([(getSubjectFilenames, calc_residuals, [('motion_param', 'motion_file')])])

    #         wf.connect([(getSubjectFilenames, globalSignalRemoval, [('brain','in_file')] )])
    #         wf.connect([(getSubjectFilenames, globalSignalRemoval, [('mask','mask_file')])])

            wf.connect([(getSubjectFilenames, bandpass, [('brain','in_file')])])
            wf.connect([(getSubjectFilenames, bandpass, [('tr','tr')])])

            wf.connect([( bandpass, pearcoff, [('out_file','in_file')])])
            wf.connect([( getSubjectFilenames, pearcoff, [('atlas','atlas_file')])])
            wf.connect([( getSubjectFilenames, pearcoff, [('mask','mask_file')])])


            # ---------------------------------------------------------------------------------------
            wf.connect([(pearcoff, func2std_xform, [('fc_map_brain_file','in_file')])])
            wf.connect([(getSubjectFilenames, func2std_xform, [('func2std_mat','in_matrix_file')])])
            wf.connect([(getSubjectFilenames, func2std_xform, [('MNI3mm_path','reference')])])

    #         -- send out file to save file list and then save the outputs



            folder_name = 'pearcoff_' + combination + '.@fc_map_brain_file'



            wf.connect([(func2std_xform,  save_file_list, [('out_file','in_fc_map_brain_file')])])
            # --------------------------------------------------------------------------------------------

            wf.connect([(save_file_list,  dataSink, [('out_fc_map_brain_file',folder_name)])])

            wf.write_graph(graph2use='flat', format='png')
            wf.run('MultiProc', plugin_args={'n_procs': num_proc})

    else:

            # wf.connect([( getSubjectFilenames, pearcoff, [('brain','in_file')])])

            # ---------------------------------------------------------------------------------------
            wf.connect([(getSubjectFilenames, highpass, [('brain','in_file')])])
            wf.connect([(getSubjectFilenames, highpass, [('tr','tr')])])

            wf.connect([( highpass, pearcoff, [('out_file','in_file')])])
            wf.connect([( getSubjectFilenames, pearcoff, [('atlas','atlas_file')])])
            wf.connect([( getSubjectFilenames, pearcoff, [('mask','mask_file')])])


            wf.connect([(pearcoff, func2std_xform, [('fc_map_brain_file','in_file')])])
            wf.connect([(getSubjectFilenames, func2std_xform, [('func2std_mat','in_matrix_file')])])
            wf.connect([(getSubjectFilenames, func2std_xform, [('MNI3mm_path','reference')])])

    #         -- send out file to save file list and then save the outputs



            folder_name = 'pearcoff_' + combination + '.@fc_map_brain_file'



            wf.connect([(func2std_xform,  save_file_list, [('out_file','in_fc_map_brain_file')])])
            # --------------------------------------------------------------------------------------------


            wf.connect([(save_file_list,  dataSink, [('out_fc_map_brain_file',folder_name)])])

            wf.write_graph(graph2use='flat', format='png')
            wf.run('MultiProc', plugin_args={'n_procs': num_proc})
