
# coding: utf-8

# ### Preprocessing Pipeline
# 1. Create a BIDSDataGrabber Node to read data files
# 2. Create a IdentityInterface - infosource Node to iterate over multiple Subjects
# 3. Create following Nodes for preprocessing
#     - [x] Exclude 4 volumes from the functional scan
#     - [x] slice time correction
#     - [x] motion correction and saving the motion parameters
#     - [x] Registration of functional data to anatomical and anatomical to standard space to create
#           transformation matrices.
#     - [x] Registering the atlas to the functional space



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


# import logging
#
# logger = logging.getLogger(__name__)
# logger.setLevel(logging.DEBUG)
#
# # create a file handler
# handler = logging.FileHandler('progress.log')
#
# # add the handlers to the logger
# logger.addHandler(handler)





# Paths

# path_cwd = os.getcwd()
# path_split_list = path_cwd.split('/')
# s = path_split_list[0:-1] # for getting to the parent dir of pwd
# s = opj('/',*s) # *s converts list to path, # very important to add '/' in the begining so it is read as directory later
#
#
# # json_path = opj(data_directory,'task-rest_bold.json')
#
# json_path = 'scripts/json/paths.json'
# with open(json_path, 'rt') as fp:
#     task_info = json.load(fp)
#
#
#
# # In[851]:
#
#
# base_directory = opj(s,task_info["base_directory_for_results"])
# parent_wf_directory = task_info["parent_wf_directory"]
# motion_correction_bet_directory = task_info["motion_correction_bet_directory"]
# coreg_reg_directory = task_info["coreg_reg_directory"]
# atlas_resize_reg_directory = task_info["atlas_resize_reg_directory"]
# data_directory = opj(s,task_info["data_directory"])
# datasink_name = task_info["datasink_name"]
#
# atlasPath = opj(s,task_info["atlas_path"])


#
# layout = BIDSLayout(data_directory)
#
# # number_of_subjects = 4 # Number of subjects you wish to preprocess
#
# subject_list = (layout.get_subjects())[0:number_of_subjects]


def main(paths, options_binary_string, ANAT , num_proc = 7):

    json_path=paths[0]
    base_directory=paths[1]
    motion_correction_bet_directory=paths[2]
    parent_wf_directory=paths[3]
    # functional_connectivity_directory=paths[4]
    coreg_reg_directory=paths[5]
    atlas_resize_reg_directory=paths[6]
    subject_list = paths[7]
    datasink_name=paths[8]
    # fc_datasink_name=paths[9]
    atlasPath=paths[10]
    # brain_path=paths[11]
    # mask_path=paths[12]
    # atlas_path=paths[13]
    # tr_path=paths[14]
    # motion_params_path=paths[15]
    # func2std_mat_path=paths[16]
    # MNI3mm_path=paths[17]
    # demographics_file_path = paths[18]
    # phenotype_file_path = paths[19]
    data_directory = paths[20]



    number_of_subjects = len(subject_list)
    print("Working with ",number_of_subjects," subjects.")

    # Create our own custom function - BIDSDataGrabber using a Function Interface.

    # In[858]:


    def get_nifti_filenames(subject_id,data_dir):
    #     Remember that all the necesary imports need to be INSIDE the function for the Function Interface to work!
        from bids.grabbids import BIDSLayout

        layout = BIDSLayout(data_dir)
        run = 1

        anat_file_path = [f.filename for f in layout.get(subject=subject_id, type='T1w', extensions=['nii', 'nii.gz'])]
        func_file_path = [f.filename for f in layout.get(subject=subject_id, type='bold', run=run, extensions=['nii', 'nii.gz'])]

        if len(anat_file_path) == 0:
            return None, func_file_path[0] # No Anatomical files present
        return anat_file_path[0],func_file_path[0]


    BIDSDataGrabber = Node(Function(function=get_nifti_filenames, input_names=['subject_id','data_dir'],
                                    output_names=['anat_file_path','func_file_path']), name='BIDSDataGrabber')
    # BIDSDataGrabber.iterables = [('subject_id',subject_list)]
    BIDSDataGrabber.inputs.data_dir = data_directory


    # ## Return TR

    def get_TR(in_file):
        from bids.grabbids import BIDSLayout

        data_directory = '/home1/varunk/data/ABIDE1/RawDataBIDs'
        layout = BIDSLayout(data_directory)
        metadata = layout.get_metadata(path=in_file)
        TR  = metadata['RepetitionTime']
        return TR


    # ---------------- Added new Node to return TR and other slice timing correction params-------------------------------
    def _getMetadata(in_file):
        from bids.grabbids import BIDSLayout
        import logging

        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)

        # create a file handler
        handler = logging.FileHandler('progress.log')

        # add the handlers to the logger
        logger.addHandler(handler)

        interleaved = True
        index_dir = False
        data_directory = '/home1/varunk/data/ABIDE1/RawDataBIDs'
        layout = BIDSLayout(data_directory)
        metadata = layout.get_metadata(path=in_file)
        print(metadata)

        logger.info('Extracting Meta Data of file: %s',in_file)
        try: tr  = metadata['RepetitionTime']
        except KeyError:
            print('Key RepetitionTime not found in task-rest_bold.json so using a default of 2.0 ')
            tr = 2
            logger.error('Key RepetitionTime not found in task-rest_bold.json for file %s so using a default of 2.0 ', in_file)

        try: slice_order = metadata['SliceAcquisitionOrder']
        except KeyError:
            print('Key SliceAcquisitionOrder not found in task-rest_bold.json so using a default of interleaved ascending ')
            logger.error('Key SliceAcquisitionOrder not found in task-rest_bold.json for file %s so using a default of interleaved ascending', in_file)
            return tr, index_dir, interleaved


        if slice_order.split(' ')[0] == 'Sequential':
            interleaved =  False
        if slice_order.split(' ')[1] == 'Descending':
            index_dir = True

        return tr, index_dir, interleaved


    getMetadata = Node(Function(function=_getMetadata, input_names=['in_file'],
                                    output_names=['tr','index_dir','interleaved']), name='getMetadata')

    # ### Skipping 4 starting scans
    # Extract ROI for skipping first 4 scans of the functional data
    # > **Arguments:**
    # t_min: (corresponds to time dimension) Denotes the starting time of the inclusion
    # t_size: Denotes the number of scans to include
    #
    # The logic behind skipping 4 initial scans is to take scans after the subject has stabalized in the scanner.

    # In[863]:


    # ExtractROI - skip dummy scans
    extract = Node(ExtractROI(t_min=4, t_size=-1),
                   output_type='NIFTI',
                   name="extract")


    # ### Slice time correction
    # Created a Node that does slice time correction
    # > **Arguments**:
    # index_dir=False -> Slices were taken bottom to top i.e. in ascending order
    # interleaved=True means odd slices were acquired first and then even slices [or vice versa(Not sure)]

    slicetimer = Node(SliceTimer(
                                 output_type='NIFTI'
                                 ),
                      name="slicetimer")


    # ### Motion Correction
    # Motion correction is done using fsl's mcflirt. It alligns all the volumes of a functional scan to each other

    # MCFLIRT - motion correction
    mcflirt = Node(MCFLIRT( mean_vol=True,
                           save_plots=True,
                           output_type='NIFTI'),
                   name="mcflirt")


    #  Just a dummy node to transfer the output of Mcflirt to the next workflow. Needed if we didnt want to use the Mcflirt
    from_mcflirt = Node(IdentityInterface(fields=['in_file']),
                      name="from_mcflirt")


    # ### Skull striping
    # I used fsl's BET

    # In[868]:


    skullStrip = Node(BET(mask=False, frac=0.3, robust=True ),name='skullStrip') #


    # *Note*: Do not include special characters in ```name``` field above coz then  wf.writegraph will cause issues

    # ## Resample
    # I needed to resample the anatomical file from 1mm to 3mm. Because registering a 1mm file was taking a huge amount of time.
    #

    # In[872]:


    # Resample - resample anatomy to 3x3x3 voxel resolution
    resample_mni = Node(Resample(voxel_size=(3, 3, 3), resample_mode='Cu', # cubic interpolation
                             outputtype='NIFTI'),
                    name="resample_mni")

    resample_anat = Node(Resample(voxel_size=(3, 3, 3), resample_mode='Cu', # cubic interpolation
                             outputtype='NIFTI'),
                    name="resample_anat")


    # In[873]:


    resample_atlas = Node(Resample(voxel_size=(3, 3, 3), resample_mode='NN', # cubic interpolation
                             outputtype='NIFTI'),
                    name="resample_atlas")

    resample_atlas.inputs.in_file = atlasPath



    # # Matrix operations
    # ### For concatenating the transformation matrices


    concat_xform = Node(ConvertXFM(concat_xfm=True),name='concat_xform')


    # Node to calculate the inverse of func2std matrix
    inv_mat = Node(ConvertXFM(invert_xfm=True), name='inv_mat')

    # ## Extracting the mean brain

    meanfunc = Node(interface=ImageMaths(op_string='-Tmean',
                                                suffix='_mean'),
                       name='meanfunc')


    meanfuncmask = Node(interface=BET(mask=True,
                                             no_output=True,
                                             frac=0.3),
                           name='meanfuncmask')

    # ## Apply Mask


    # Does BET (masking) on the whole func scan [Not using this, creates bug for join node]
    maskfunc = Node(interface=ImageMaths(suffix='_bet',
                                                   op_string='-mas'),
                          name='maskfunc')

    # Does BET (masking) on the mean func scan
    maskfunc4mean = Node(interface=ImageMaths(suffix='_bet',
                                                   op_string='-mas'),
                          name='maskfunc4mean')

    # ## Datasink
    # I needed to define the structure of what files are saved and where.

    # Create DataSink object
    dataSink = Node(DataSink(), name='datasink')

    # Name of the output folder
    dataSink.inputs.base_directory = opj(base_directory,datasink_name)


    # Define substitution strings so that the data is similar to BIDS
    substitutions = [('_subject_id_', 'sub-'),
                     ('_resample_brain_flirt.nii_brain', ''),
                     ('_roi_st_mcf_flirt.nii_brain_flirt', ''),
                     ('task-rest_run-1_bold_roi_st_mcf.nii','motion_params'),
                     ('T1w_resample_brain_flirt_sub-0050002_task-rest_run-1_bold_roi_st_mcf_mean_bet_flirt','fun2std')
                    ]

    # Feed the substitution strings to the DataSink node
    dataSink.inputs.substitutions = substitutions


    # ### Apply Mask to functional data
    # Mean file of the motion corrected functional scan is sent to
    # skullStrip to get just the brain and the mask_image.
    # Mask_image is just a binary file (containing 1 where brain is present and 0 where it isn't).
    # After getting the mask_image form skullStrip, apply that mask to aligned
    # functional image to extract its brain and remove the skull

    # In[889]:


    # Function
    # in_file: The file on which you want to apply mask
    # in_file2 = mask_file:  The mask you want to use. Make sure that mask_file has same size as in_file
    # out_file : Result of applying mask in in_file -> Gives the path of the output file

    def applyMask_func(in_file, in_file2):
        import numpy as np
        import nibabel as nib
        import os
        from os.path import join as opj

        # convert from unicode to string : u'/tmp/tmp8daO2Q/..' -> '/tmp/tmp8daO2Q/..' i.e. removes the prefix 'u'
        mask_file = in_file2

        brain_data = nib.load(in_file)
        mask_data = nib.load(mask_file)

        brain = brain_data.get_data().astype('float32')
        mask = mask_data.get_data()

        # applying mask by multiplying elementwise to the binary mask

        if len(brain.shape) == 3: # Anat file
            brain = np.multiply(brain,mask)
        elif len(brain.shape) > 3: # Functional File
            for t in range(brain.shape[-1]):
                brain[:,:,:,t] = np.multiply(brain[:,:,:,t],mask)
        else:
            pass

        # Saving the brain file

        path = os.getcwd()


        in_file_split_list = in_file.split('/')
        in_file_name = in_file_split_list[-1]

        out_file = in_file_name + '_brain.nii.gz' # changing name
        brain_with_header = nib.Nifti1Image(brain, affine=brain_data.affine,header = brain_data.header)
        nib.save(brain_with_header,out_file)

        out_file = opj(path,out_file)
        out_file2 = in_file2

        return out_file, out_file2



    # #### Things learnt:
    # 1. I found out that whenever a node is being executed, it becomes the current directory and whatever file you create now, will be stored here.
    # 2. #from IPython.core.debugger import Tracer; Tracer()()    # Debugger doesnt work in nipype

    # Wrap the above function inside a Node

    # In[890]:


    applyMask = Node(Function(function=applyMask_func, input_names=['in_file','in_file2'],
                                    output_names=['out_file','out_file2']), name='applyMask')


    # ### Some nodes needed for Co-registration and Normalization


    # Node for getting the xformation matrix
    func2anat_reg = Node(FLIRT(output_type='NIFTI'), name="func2anat_reg")

    # Node for applying xformation matrix to functional data
    func2std_xform = Node(FLIRT(output_type='NIFTI',
                             apply_xfm=True), name="func2std_xform")

    # Node for applying xformation matrix to functional data
    std2func_xform = Node(FLIRT(output_type='NIFTI',
                             apply_xfm=True, interp='nearestneighbour'), name="std2func_xform")


    # Node for Normalizing/Standardizing the anatomical and getting the xformation matrix
    anat2std_reg = Node(FLIRT(output_type='NIFTI'), name="anat2std_reg")




    # I wanted to use the MNI file as input to the workflow so I created an Identity
    # Node that reads the MNI file path and outputs the same MNI file path.
    # Then I connected this node to whereever it was needed.



    MNI152_2mm = Node(IdentityInterface(fields=['standard_file','mask_file']),
                      name="MNI152_2mm")
    # Set the mask_file and standard_file input in the Node. This setting sets the input mask_file permanently.
    MNI152_2mm.inputs.mask_file = os.path.expandvars('$FSLDIR/data/standard/MNI152_T1_2mm_brain_mask.nii.gz')

    MNI152_2mm.inputs.standard_file = os.path.expandvars('$FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz')
    # MNI152_2mm.inputs.mask_file = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain_mask.nii.gz'
    # MNI152_2mm.inputs.standard_file = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain.nii.gz'



    # ## Band Pass Filtering
    # Let's do a band pass filtering on the data using the code from https://neurostars.org/t/bandpass-filtering-different-outputs-from-fsl-and-nipype-custom-function/824/2



    ### AFNI

    bandpass = Node(afni.Bandpass(highpass=0.008, lowpass=0.08,
                             despike=False, no_detrend=True, notrans=True,
                             outputtype='NIFTI_GZ'),name='bandpass')

    # ### Following is a Join Node that collects the preprocessed file paths and saves them in a file

    # In[902]:


    def save_file_list_function_in_brain(in_brain):
        import numpy as np
        import os
        from os.path import join as opj


        file_list = np.asarray(in_brain)
        print('######################## File List ######################: \n',file_list)

        np.save('brain_file_list',file_list)
        file_name = 'brain_file_list.npy'
        out_brain = opj(os.getcwd(),file_name) # path
        return out_brain

    def save_file_list_function_in_mask(in_mask):
          import numpy as np
          import os
          from os.path import join as opj

          file_list2 = np.asarray(in_mask)
          print('######################## File List ######################: \n',file_list2)

          np.save('mask_file_list',file_list2)
          file_name2 = 'mask_file_list.npy'
          out_mask = opj(os.getcwd(),file_name2) # path
          return out_mask

    def save_file_list_function_in_motion_params(in_motion_params):
          import numpy as np
          import os
          from os.path import join as opj

          file_list3 = np.asarray(in_motion_params)
          print('######################## File List ######################: \n',file_list3)

          np.save('motion_params_file_list',file_list3)
          file_name3 = 'motion_params_file_list.npy'
          out_motion_params = opj(os.getcwd(),file_name3) # path
          return out_motion_params

    def save_file_list_function_in_motion_outliers(in_motion_outliers):
          import numpy as np
          import os
          from os.path import join as opj

          file_list4 = np.asarray(in_motion_outliers)
          print('######################## File List ######################: \n',file_list4)

          np.save('motion_outliers_file_list',file_list4)
          file_name4 = 'motion_outliers_file_list.npy'
          out_motion_outliers = opj(os.getcwd(),file_name4) # path
          return out_motion_outliers

    def save_file_list_function_in_joint_xformation_matrix(in_joint_xformation_matrix):
          import numpy as np
          import os
          from os.path import join as opj

          file_list5 = np.asarray(in_joint_xformation_matrix)
          print('######################## File List ######################: \n',file_list5)

          np.save('joint_xformation_matrix_file_list',file_list5)
          file_name5 = 'joint_xformation_matrix_file_list.npy'
          out_joint_xformation_matrix = opj(os.getcwd(),file_name5) # path
          return out_joint_xformation_matrix

    def save_file_list_function_in_tr(in_tr):
          import numpy as np
          import os
          from os.path import join as opj

          tr_list = np.asarray(in_tr)
          print('######################## TR List ######################: \n',tr_list)

          np.save('tr_list',tr_list)
          file_name6 = 'tr_list.npy'
          out_tr = opj(os.getcwd(),file_name6) # path
          return out_tr

    def save_file_list_function_in_atlas(in_atlas):
          import numpy as np
          import os
          from os.path import join as opj

          file_list7 = np.asarray(in_atlas)
          print('######################## File List ######################: \n',file_list7)

          np.save('atlas_file_list',file_list7)
          file_name7 = 'atlas_file_list.npy'
          out_atlas = opj(os.getcwd(),file_name7) # path
          return out_atlas


    save_file_list_in_brain = JoinNode(Function(function=save_file_list_function_in_brain, input_names=['in_brain'],
                   output_names=['out_brain']),
                   joinsource="infosource",
                   joinfield=['in_brain'],
                   name="save_file_list_in_brain")

    save_file_list_in_mask = JoinNode(Function(function=save_file_list_function_in_mask, input_names=['in_mask'],
                   output_names=['out_mask']),
                   joinsource="infosource",
                   joinfield=['in_mask'],
                   name="save_file_list_in_mask")


    save_file_list_in_motion_outliers = JoinNode(Function(function=save_file_list_function_in_motion_outliers, input_names=['in_motion_outliers'],
                   output_names=['out_motion_outliers']),
                   joinsource="infosource",
                   joinfield=['in_motion_outliers'],
                   name="save_file_list_in_motion_outliers")


    save_file_list_in_motion_params = JoinNode(Function(function=save_file_list_function_in_motion_params, input_names=['in_motion_params'],
                   output_names=['out_motion_params']),
                   joinsource="infosource",
                   joinfield=['in_motion_params'],
                   name="save_file_list_in_motion_params")



    save_file_list_in_joint_xformation_matrix = JoinNode(Function(function=save_file_list_function_in_joint_xformation_matrix, input_names=['in_joint_xformation_matrix'],
                   output_names=['out_joint_xformation_matrix']),
                   joinsource="infosource",
                   joinfield=['in_joint_xformation_matrix'],
                   name="save_file_list_in_joint_xformation_matrix")

    save_file_list_in_tr = JoinNode(Function(function=save_file_list_function_in_tr, input_names=['in_tr'],
               output_names=['out_tr']),
               joinsource="infosource",
               joinfield=['in_tr'],
               name="save_file_list_in_tr")



    save_file_list_in_atlas = JoinNode(Function(function=save_file_list_function_in_atlas, input_names=['in_atlas'],
               output_names=['out_atlas']),
               joinsource="infosource",
               joinfield=['in_atlas'],
               name="save_file_list_in_atlas")


    # save_file_list = JoinNode(Function(function=save_file_list_function, input_names=['in_brain', 'in_mask', 'in_motion_params','in_motion_outliers','in_joint_xformation_matrix', 'in_tr', 'in_atlas'],
    #                output_names=['out_brain','out_mask','out_motion_params','out_motion_outliers','out_joint_xformation_matrix','out_tr', 'out_atlas']),
    #                joinsource="infosource",
    #                joinfield=['in_brain', 'in_mask', 'in_motion_params','in_motion_outliers','in_joint_xformation_matrix','in_tr', 'in_atlas'],
    #                name="save_file_list")








    # def save_file_list_function(in_brain, in_mask, in_motion_params, in_motion_outliers, in_joint_xformation_matrix, in_tr, in_atlas):
    #     # Imports
    #     import numpy as np
    #     import os
    #     from os.path import join as opj
    #
    #
    #     file_list = np.asarray(in_brain)
    #     print('######################## File List ######################: \n',file_list)
    #
    #     np.save('brain_file_list',file_list)
    #     file_name = 'brain_file_list.npy'
    #     out_brain = opj(os.getcwd(),file_name) # path
    #
    #
    #     file_list2 = np.asarray(in_mask)
    #     print('######################## File List ######################: \n',file_list2)
    #
    #     np.save('mask_file_list',file_list2)
    #     file_name2 = 'mask_file_list.npy'
    #     out_mask = opj(os.getcwd(),file_name2) # path
    #
    #
    #     file_list3 = np.asarray(in_motion_params)
    #     print('######################## File List ######################: \n',file_list3)
    #
    #     np.save('motion_params_file_list',file_list3)
    #     file_name3 = 'motion_params_file_list.npy'
    #     out_motion_params = opj(os.getcwd(),file_name3) # path
    #
    #
    #     file_list4 = np.asarray(in_motion_outliers)
    #     print('######################## File List ######################: \n',file_list4)
    #
    #     np.save('motion_outliers_file_list',file_list4)
    #     file_name4 = 'motion_outliers_file_list.npy'
    #     out_motion_outliers = opj(os.getcwd(),file_name4) # path
    #
    #
    #     file_list5 = np.asarray(in_joint_xformation_matrix)
    #     print('######################## File List ######################: \n',file_list5)
    #
    #     np.save('joint_xformation_matrix_file_list',file_list5)
    #     file_name5 = 'joint_xformation_matrix_file_list.npy'
    #     out_joint_xformation_matrix = opj(os.getcwd(),file_name5) # path
    #
    #     tr_list = np.asarray(in_tr)
    #     print('######################## TR List ######################: \n',tr_list)
    #
    #     np.save('tr_list',tr_list)
    #     file_name6 = 'tr_list.npy'
    #     out_tr = opj(os.getcwd(),file_name6) # path
    #
    #
    #     file_list7 = np.asarray(in_atlas)
    #     print('######################## File List ######################: \n',file_list7)
    #
    #     np.save('atlas_file_list',file_list7)
    #     file_name7 = 'atlas_file_list.npy'
    #     out_atlas = opj(os.getcwd(),file_name7) # path
    #
    #
    #
    #
    #     return out_brain, out_mask, out_motion_params, out_motion_outliers, out_joint_xformation_matrix, out_tr , out_atlas
    #
    #
    #
    # save_file_list = JoinNode(Function(function=save_file_list_function, input_names=['in_brain', 'in_mask', 'in_motion_params','in_motion_outliers','in_joint_xformation_matrix', 'in_tr', 'in_atlas'],
    #                  output_names=['out_brain','out_mask','out_motion_params','out_motion_outliers','out_joint_xformation_matrix','out_tr', 'out_atlas']),
    #                  joinsource="infosource",
    #                  joinfield=['in_brain', 'in_mask', 'in_motion_params','in_motion_outliers','in_joint_xformation_matrix','in_tr', 'in_atlas'],
    #                  name="save_file_list")


    # ### Motion outliers


    motionOutliers = Node(MotionOutliers(no_motion_correction=False,metric='fd', out_metric_plot = 'fd_plot.png',
                                         out_metric_values='fd_raw.txt'),name='motionOutliers')

    # ## Workflow for atlas registration  from std to functional

    wf_atlas_resize_reg = Workflow(name=atlas_resize_reg_directory)

    wf_atlas_resize_reg.connect([

                # Apply the inverse matrix to the 3mm Atlas to transform it to func space

                (maskfunc4mean, std2func_xform, [(('out_file','reference'))]),

                (resample_atlas, std2func_xform, [('out_file','in_file')] ),

                # Now, applying the inverse matrix

                (inv_mat, std2func_xform, [('out_file','in_matrix_file')]), # output: Atlas in func space

                (std2func_xform, save_file_list_in_atlas, [('out_file','in_atlas')]),

                # ---------------------------Save the required files --------------------------------------------

                (save_file_list_in_motion_params, dataSink, [('out_motion_params','motion_params_paths.@out_motion_params')]),
                (save_file_list_in_motion_outliers, dataSink, [('out_motion_outliers','motion_outliers_paths.@out_motion_outliers')]),
                (save_file_list_in_brain, dataSink, [('out_brain','preprocessed_brain_paths.@out_brain')]),
                (save_file_list_in_mask, dataSink, [('out_mask','preprocessed_mask_paths.@out_mask')]),

                (save_file_list_in_joint_xformation_matrix, dataSink, [('out_joint_xformation_matrix',
                                             'joint_xformation_matrix_paths.@out_joint_xformation_matrix')]),

                (save_file_list_in_tr, dataSink, [('out_tr','tr_paths.@out_tr')]),

                (save_file_list_in_atlas, dataSink, [('out_atlas','atlas_paths.@out_atlas')])



    ])


    # In[909]:


    wf_coreg_reg = Workflow(name=coreg_reg_directory)
    # wf_coreg_reg.base_dir = base_directory
    # Dir where all the outputs will be stored(inside coregistrationPipeline folder).



    if ANAT == 1:
        wf_coreg_reg.connect(BIDSDataGrabber,'anat_file_path',skullStrip,'in_file') # Resampled the anat file to 3mm

        wf_coreg_reg.connect(skullStrip,'out_file', resample_anat,'in_file')

        wf_coreg_reg.connect(resample_anat,'out_file', func2anat_reg,'reference') # Make the resampled file as reference in func2anat_reg

        # Sec 1. The above 3 steps registers the mean image to resampled anat image and
        # calculates the xformation matrix .. I hope the xformation matrix will be saved

        wf_coreg_reg.connect(MNI152_2mm, 'standard_file', resample_mni,'in_file')

        wf_coreg_reg.connect(resample_mni, 'out_file', anat2std_reg,'reference')

        wf_coreg_reg.connect(resample_anat, 'out_file', anat2std_reg, 'in_file')


        # Calculates the Xformationmatrix from anat3mm to MNI 3mm

        # We can get those matrices by refering to func2anat_reg.outputs.out_matrix_file and similarly for anat2std_reg





        wf_coreg_reg.connect(func2anat_reg, 'out_matrix_file', concat_xform,'in_file')

        wf_coreg_reg.connect(anat2std_reg, 'out_matrix_file', concat_xform,'in_file2')

        wf_coreg_reg.connect(concat_xform, 'out_file', dataSink, 'tranformation_matrix_fun2std.@out_file')

        wf_coreg_reg.connect(concat_xform, 'out_file', save_file_list_in_joint_xformation_matrix, 'in_joint_xformation_matrix')

        # Now inverse the func2std MAT to std2func
        wf_coreg_reg.connect(concat_xform, 'out_file', wf_atlas_resize_reg,'inv_mat.in_file')
# ------------------------------------------------------------------------------------------------------------------------------

    # Registration of Functional to MNI 3mm space w/o using anatomical
    if ANAT == 0:
        print('Not using Anatomical high resoulution files')
        wf_coreg_reg.connect(MNI152_2mm, 'standard_file', resample_mni,'in_file')
        wf_coreg_reg.connect(resample_mni, 'out_file',func2anat_reg,'reference') # Make the resampled file as reference in func2anat_reg

        wf_coreg_reg.connect(func2anat_reg, 'out_matrix_file', dataSink, 'tranformation_matrix_fun2std.@out_file')

        wf_coreg_reg.connect(func2anat_reg, 'out_matrix_file', save_file_list_in_joint_xformation_matrix, 'in_joint_xformation_matrix')

        # Now inverse the func2std MAT to std2func
        wf_coreg_reg.connect(func2anat_reg, 'out_matrix_file', wf_atlas_resize_reg,'inv_mat.in_file')









    # ## Co-Registration, Normalization and Bandpass Workflow
    # 1. Co-registration means alligning the func to anat
    # 2. Normalization means aligning func/anat to standard
    # 3. Applied band pass filtering in range - highpass=0.008, lowpass=0.08

    # In[910]:


    wf_motion_correction_bet = Workflow(name=motion_correction_bet_directory)
    # wf_motion_correction_bet.base_dir = base_directory

    wf_motion_correction_bet.connect([


                (from_mcflirt, meanfunc, [('in_file','in_file')]),

                (meanfunc, meanfuncmask, [('out_file','in_file')]),




                (from_mcflirt,applyMask , [('in_file','in_file')]), # 1

                (meanfuncmask, applyMask, [('mask_file','in_file2')]), # 2 output: 1&2,  BET on coregistered fmri scan

                (meanfunc, maskfunc4mean, [('out_file', 'in_file')]), # 3

                (meanfuncmask, maskfunc4mean, [('mask_file','in_file2')]), # 4 output: 3&4, BET on mean func scan


                (applyMask, save_file_list_in_brain, [('out_file', 'in_brain')]),
                (applyMask, save_file_list_in_mask, [('out_file2', 'in_mask')]),



                (maskfunc4mean, wf_coreg_reg, [('out_file','func2anat_reg.in_file')])



    ])



    infosource = Node(IdentityInterface(fields=['subject_id']),
                      name="infosource")

    infosource.iterables = [('subject_id',subject_list)]



    # Create the workflow

    wf = Workflow(name=parent_wf_directory)
    # base_dir = opj(s,'result')
    wf.base_dir = base_directory # Dir where all the outputs will be stored(inside BETFlow folder).

    # wf.connect([      (infosource, BIDSDataGrabber, [('subject_id','subject_id')]),
    #                   (BIDSDataGrabber, extract, [('func_file_path','in_file')]),
    #
    #                   (BIDSDataGrabber,getMetadata, [('func_file_path','in_file')]),
    #
    #                   (getMetadata,slicetimer, [('tr','time_repetition')]),
    #
    #
    #                   (getMetadata,slicetimer, [('index_dir','index_dir')]),
    #
    #                   (getMetadata,slicetimer, [('interleaved','interleaved')]),
    #
    #                   (getMetadata,save_file_list_in_tr, [('tr','in_tr')]),
    #
    #                   (extract,slicetimer,[('roi_file','in_file')]),
    #
    #                   (slicetimer, mcflirt,[('slice_time_corrected_file','in_file')])
    #                   (mcflirt,dataSink,[('par_file','motion_params.@par_file')]), # saves the motion parameters calculated before
    #
    #                   (mcflirt,save_file_list_in_motion_params,[('par_file','in_motion_params')]),
    #
    #                   (mcflirt,wf_motion_correction_bet,[('out_file','from_mcflirt.in_file')])
    #            ])
    # # Run it in parallel
    # wf.run('MultiProc', plugin_args={'n_procs': num_proc})
    #
    #
    #
    # # Visualize the detailed graph
    # # from IPython.display import Image
    # wf.write_graph(graph2use='flat', format='png', simple_form=True)



    # Options:
    # discard 4 Volumes (extract), slicetimer, mcflirt
    print('Preprocessing Options:')
    print('Skipping 4 dummy volumes - ',options_binary_string[0])
    print('Slicetiming correction - ',options_binary_string[1])
    print('Finding Motion Outliers - ',options_binary_string[2])
    print('Doing Motion Correction - ',options_binary_string[3])

    # ANAT = 0
    nodes = [extract, slicetimer,motionOutliers, mcflirt]
    wf.connect(infosource,'subject_id', BIDSDataGrabber,'subject_id')
    wf.connect(BIDSDataGrabber, 'func_file_path', getMetadata, 'in_file')
    wf.connect(getMetadata, 'tr', save_file_list_in_tr,'in_tr')

    old_node = BIDSDataGrabber
    old_node_output = 'func_file_path'

    for idx, include  in enumerate(options_binary_string):

        if old_node == extract :
            old_node_output = 'roi_file'
        elif old_node == slicetimer:
            old_node_output = 'slice_time_corrected_file'
        # elif old_node == mcflirt:



            # old_node_output = 'out_file'



        if int(include):
            new_node = nodes[idx]



            if new_node == slicetimer:
                wf.connect(getMetadata,'tr',slicetimer,'time_repetition')
                wf.connect(getMetadata,'index_dir',slicetimer, 'index_dir')
                wf.connect(getMetadata,'interleaved',slicetimer,'interleaved')
                new_node_input = 'in_file'
            elif new_node == extract:
                new_node_input = 'in_file'
            elif new_node == mcflirt:
                new_node_input = 'in_file'
                wf.connect(mcflirt,'par_file',dataSink,'motion_params.@par_file') # saves the motion parameters calculated before

                wf.connect(mcflirt,'par_file',save_file_list_in_motion_params,'in_motion_params')

                wf.connect(mcflirt, 'out_file', wf_motion_correction_bet,'from_mcflirt.in_file')


            elif new_node == motionOutliers:

                wf.connect(meanfuncmask, 'mask_file', motionOutliers,'mask')

                wf.connect(motionOutliers, 'out_file', dataSink,'motionOutliers.@out_file')

                wf.connect(motionOutliers, 'out_metric_plot', dataSink,'motionOutliers.@out_metric_plot')

                wf.connect(motionOutliers, 'out_metric_values', dataSink,'motionOutliers.@out_metric_values')

                wf.connect(motionOutliers, 'out_file', save_file_list_in_motion_outliers,'in_motion_outliers')

                new_node_input = 'in_file'

                wf.connect(old_node, old_node_output, new_node, new_node_input)


                continue



            wf.connect(old_node, old_node_output, new_node, new_node_input)

            old_node = new_node

        else:
            if idx == 3:
                # new_node = from_mcflirt
                # new_node_input = 'from_mcflirt.in_file'

                wf.connect(old_node, old_node_output, wf_motion_correction_bet,'from_mcflirt.in_file')

                # old_node = new_node




    TEMP_DIR_FOR_STORAGE = opj(base_directory,'crash_files')
    wf.config = {"execution": {"crashdump_dir": TEMP_DIR_FOR_STORAGE}}

    # Visualize the detailed graph
    # from IPython.display import Image

    wf.write_graph(graph2use='flat', format='png', simple_form=True)

    # Run it in parallel
    wf.run('MultiProc', plugin_args={'n_procs': num_proc})
