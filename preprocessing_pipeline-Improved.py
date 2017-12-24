
# coding: utf-8

# ### Preprocessing Pipeline
# 1. Create a BIDSDataGrabber Node to read data files
# 2. Create a IdentityInterface Node to iterate over multiple Subjects
# 3. Create following Nodes for preprocessing: (Based on [Nan-kuei Chen's resting state analysis pipeline:](https://wiki.biac.duke.edu/biac:analysis:resting_pipeline)
#     - [-] convert data to nii in LAS orientation (Skip if NYU is already in [LAS Orientation](http://www.grahamwideman.com/gw/brain/orientation/orientterms.htm))
#     - [x] Exclude 4 volumes from the functional scan
#     - [x] slice time correction
#     - [x] motion correction, {[then regress out motion parameter] - This will be done later}
#     - [x] Skull stripping and mask generation using mean of functional scan got using mcflirt
#     - [x] Apply mask to Functional image
#     - [x] Co-Registration with Anatomical Image
#     - [x] normalize functional data
#     - [-] regress out WM/CSF - Not doing coz of the debate that WM also has some activations
#     - [x] bandpass filter
#
# 4. Embed them into a workflow
# 5. Do the Preprocessing of 4 subjects

# In[846]:


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


# ## Define Paths
# Let's set the directory names:
# 1. **base_directory** : The directory where all the output of my program will be saved
# 2. I have created 2 workflows, one onside another:
#     3. **parent_wf_directory**: The name of the folder where the top level workflow's output is saved
#     4. **child_wf_directory**: The name of the folder where the Second level workfolw's output is saved
# 5. **data_directory**: Directory where the BIDS data is stored.
#

# In[847]:


# get_ipython().system('pwd')


# ## Adding module to read the parameters and paths from json file

# In[848]:





# In[849]:


# Paths

path_cwd = os.getcwd()
path_split_list = path_cwd.split('/')
s = path_split_list[0:-1] # for getting to the parent dir of pwd
s = opj('/',*s) # *s converts list to path, # very important to add '/' in the begining so it is read as directory later



# In[850]:



# json_path = opj(data_directory,'task-rest_bold.json')

json_path = 'scripts/json/paths.json'
with open(json_path, 'rt') as fp:
    task_info = json.load(fp)



# In[851]:


# base_directory = opj(s,'result')
# parent_wf_directory = 'preprocessPipeline_ABIDE2_GU1_withfloat'
# child_wf_directory = 'coregistrationPipeline'

# data_directory = opj(s,"data/ABIDE2-BIDS/GU1")

# datasink_name = 'datasink_preprocessed_ABIDE2_GU1_withfloat'

base_directory = opj(s,task_info["base_directory_for_results"])
motion_correction_bet_directory = task_info["motion_correction_bet_directory"]
parent_wf_directory = task_info["parent_wf_directory"]
coreg_reg_directory = task_info["coreg_reg_directory"]
atlas_resize_reg_directory = task_info["atlas_resize_reg_directory"]
data_directory = opj(s,task_info["data_directory"])
datasink_name = task_info["datasink_name"]

atlasPath = opj(s,task_info["atlas_path"])


# mask_file = '/media/varun/LENOVO4/Projects/result/preprocessPipeline/coregistrationPipeline/_subject_id_0050952/skullStrip/sub-0050952_T1w_resample_brain_mask.nii.gz'
# os.chdir(path)


# In[852]:


# path_cwd


# In[853]:


layout = BIDSLayout(data_directory)

# number_of_subjects = 4 # Number of subjects you wish to preprocess
number_of_subjects = len(layout.get_subjects())
print("Working with ",number_of_subjects," subjects.")

# Checking the Data directory Structure

# In[854]:


# !tree /home/jovyan/work/preprocess/data/ABIDE-BIDS/NYU/


# In[855]:


# len(layout.get_subjects()) # working!Gives us list of all the subjects


# In[856]:


# layout.get_subjects();


# To get the metadata associated with a subject. [Takes as argumment the filename of subject ]

# Create a list of subjects

# In[857]:


subject_list = (layout.get_subjects())[0:number_of_subjects]


# Create our own custom function - BIDSDataGrabber using a Function Interface.

# In[858]:


def get_nifti_filenames(subject_id,data_dir):
#     Remember that all the necesary imports need to be INSIDE the function for the Function Interface to work!
    from bids.grabbids import BIDSLayout

    layout = BIDSLayout(data_dir)

    anat_file_path = [f.filename for f in layout.get(subject=subject_id, type='T1w', extensions=['nii', 'nii.gz'])]
    func_file_path = [f.filename for f in layout.get(subject=subject_id, type='bold', run='1', extensions=['nii', 'nii.gz'])]

    return anat_file_path[0],func_file_path[0]

# Refer to Supplementary material section One for info on arguments for layout.get()


# Wrap it inside a Node

# In[859]:


BIDSDataGrabber = Node(Function(function=get_nifti_filenames, input_names=['subject_id','data_dir'],
                                output_names=['anat_file_path','func_file_path']), name='BIDSDataGrabber')
# BIDSDataGrabber.iterables = [('subject_id',subject_list)]
BIDSDataGrabber.inputs.data_dir = data_directory

# %%bash
#
# cd /home1/varunk/Autism-Connectome-Analysis-bids-related
# In[860]:


# To test the function wrapped in the node
# os.chdir('/home1/varunk/Autism-Connectome-Analysis-bids-related')
# BIDSDataGrabber.inputs.data_dir = data_directory
# BIDSDataGrabber.inputs.subject_id = layout.get_subjects()[0] # gives the first subject's ID
# res = BIDSDataGrabber.run()

# res.outputs


# ## Return TR

# In[861]:


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

    interleaved = True
    index_dir = False
    data_directory = '/home1/varunk/data/ABIDE1/RawDataBIDs'
    layout = BIDSLayout(data_directory)
    metadata = layout.get_metadata(path=in_file)
    tr  = metadata['RepetitionTime']
    slice_order = metadata['SliceAcquisitionOrder']
    
    if slice_order.split(' ')[0] == 'Sequential':
        interleaved =  False
    if slice_order.split(' ')[1] == 'Descending': 
        index_dir = True

    return tr, index_dir, interleaved


getMetadata = Node(Function(function=_getMetadata, input_names=['in_file'],
                                output_names=['tr','index_dir','interleaved']), name='getMetadata')

# In[862]:


# type(get_TR('/home1/varunk/data/ABIDE1/RawDataBIDs/Pitt/sub-0050002/func/sub-0050002_task-rest_run-1_bold.nii.gz'))


#
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

# In[864]:


slicetimer = Node(SliceTimer(
                             output_type='NIFTI'
                             ),
                  name="slicetimer")



# index_dir=False,interleaved=True,

# In[865]:


# To test Slicetimer

# subject_id = layout.get_subjects()[0] # gives the first subject's ID
# func_file_path = [f.filename for f in layout.get(subject=subject_id, type='bold', extensions=['nii', 'nii.gz'])]
# slicetimer.inputs.in_file = func_file_path[0]
# res = slicetimer.run()
# res.outputs


# ### Motion Correction
# Motion correction is done using fsl's mcflirt. It alligns all the volumes of a functional scan to each other

# In[866]:


# MCFLIRT - motion correction
mcflirt = Node(MCFLIRT( mean_vol=True,
                       save_plots=True,
                       output_type='NIFTI'),
               name="mcflirt")

#  ref_vol = 1,


# In[867]:


# To test mcflirt

# subject_id = layout.get_subjects()[0] # gives the first subject's ID
# func_file_path = [f.filename for f in layout.get(subject=subject_id, type='bold', extensions=['nii', 'nii.gz'])]
# mcflirt.inputs.in_file = func_file_path[0]
# res_mcflirt = mcflirt.run()
# res_mcflirt.outputs


# %%bash
# cat /tmp/tmpc2wmdeci/mcflirt/sub-28741_task-rest_run-1_bold_mcf.nii.par

# ### Skull striping
# I used fsl's BET

# In[868]:


skullStrip = Node(BET(mask=False, frac=0.3, robust=True ),name='skullStrip') # 


# *Note*: Do not include special characters in ```name``` field above coz then  wf.writegraph will cause issues

# In[869]:


# BET.help(); # Useful to see what are the parameters taken by BET


# ## Atlas

# In[870]:


# Put in the path of atlas you wish to use
# atlasPath = opj(s,'atlas/Full_brain_atlas_thr0-2mm/fullbrain_atlas_thr0-2mm.nii.gz')


# In[871]:


# # Read the atlas
# atlasObject = nib.load(atlasPath)
# atlas = atlasObject.get_data()


# ## Resample
# I needed to resample the anatomical file from 1mm to 2mm. Because registering a 1mm file was taking a huge amount of time.
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


# In[874]:


# Resample.help() # To understand what all parameters Resample supports


# In[875]:


# resample.outputs


# # Matrix operations
# ### For concatenating the transformation matrices

# In[876]:


concat_xform = Node(ConvertXFM(concat_xfm=True),name='concat_xform')
# .cmdline


# ### For finding the inverse of a transformation matrix

# In[877]:


# Node to calculate the inverse of func2std matrix
inv_mat = Node(ConvertXFM(invert_xfm=True), name='inv_mat')


# In[878]:


# inv_mat.inputs


# ## Extracting the mean brain

# In[879]:


meanfunc = Node(interface=ImageMaths(op_string='-Tmean',
                                            suffix='_mean'),
                   name='meanfunc')

# preproc.connect(motion_correct, ('out_file', pickfirst), meanfunc, 'in_file')


# In[880]:


# in_file = '/home1/varunk/data/ABIDE1/RawDataBIDs/Pitt/sub-0050002/func/sub-0050002_task-rest_run-1_bold.nii.gz'
# meanfunc.inputs.in_file = in_file
# res = meanfunc.run()


# ## Creating mask using the mean brain

# In[881]:


meanfuncmask = Node(interface=BET(mask=True,
                                         no_output=True,
                                         frac=0.3),
                       name='meanfuncmask')


# In[882]:


# in_file = '/home1/varunk/data/ABIDE1/RawDataBIDs/Pitt/sub-0050002/func/sub-0050002_task-rest_run-1_bold.nii.gz'
# meanfuncmask.inputs.in_file = in_file
# res = meanfuncmask.run()


# ## Apply Mask

# In[883]:


# Does BET (masking) on the whole func scan [Not using this, creates bug for join node]
maskfunc = Node(interface=ImageMaths(suffix='_bet',
                                               op_string='-mas'),
                      name='maskfunc')

# Does BET (masking) on the mean func scan
maskfunc4mean = Node(interface=ImageMaths(suffix='_bet',
                                               op_string='-mas'),
                      name='maskfunc4mean')


# In[884]:


# in_file = '/home1/varunk/data/ABIDE1/RawDataBIDs/Pitt/sub-0050002/func/sub-0050002_task-rest_run-1_bold.nii.gz'

# in_file = '/usr/local/fsl/data/standard/MNI152_T1_2mm.nii.gz'
# in_file2 = '/usr/local/fsl/data/standard/MNI152_T1_2mm_brain_mask.nii.gz'
# maskfunc.inputs.in_file = in_file
# maskfunc.inputs.in_file2 = in_file2
# res = maskfunc.run()


# In[885]:


# res.outputs.out_file


# ## Datasink
# I needed to define the structure of what files are saved and where.

# In[886]:


# Create DataSink object
dataSink = Node(DataSink(), name='datasink')

# Name of the output folder
dataSink.inputs.base_directory = opj(base_directory,datasink_name)




# In[887]:


# base_directory


# To create the substitutions I looked the `datasink` folder where I was redirecting the output. I manually selected the part of file/folder name that I wanted to change and copied below to be substituted.
#
# **TODO:** Using datasink create a hierarchical directory structure i.e. folder in folder - to exactly match BIDS.

# In[888]:


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
# Mean file of the motion corrected functional scan is sent to skullStrip to get just the brain and the mask_image. Mask_image is just a binary file (containing 1 where brain is present and 0 where it isn't).
# After getting the mask_image form skullStrip, apply that mask to aligned functional image to extract its brain and remove the skull

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

# I observed using fslsyes that the brain is enlarged if you Normalize a  brain resampled to 2mm brain. This in turn causes the functional data to enlarge as well after normalization. So, I will apply MNI152_2mm brain mask to the  resample brain after it has been normalized.
#
# For that let's first create a Node - `anat2std_reg_masking`  that applies the  MNI152_2mm brain mask to the Output of anat2std_reg.

# In[891]:


# FLIRT.help()


# In[892]:


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




# I wanted to use the MNI file as input to the workflow so I created an Identity Node that reads the MNI file path and outputs the same MNI file path. Then I connected this node to whereever it was needed.

# In[893]:


MNI152_2mm = Node(IdentityInterface(fields=['standard_file','mask_file']),
                  name="MNI152_2mm")
# Set the mask_file and standard_file input in the Node. This setting sets the input mask_file permanently.
MNI152_2mm.inputs.mask_file = os.path.expandvars('$FSLDIR/data/standard/MNI152_T1_2mm_brain_mask.nii.gz')

MNI152_2mm.inputs.standard_file = os.path.expandvars('$FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz')
# MNI152_2mm.inputs.mask_file = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain_mask.nii.gz'
# MNI152_2mm.inputs.standard_file = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain.nii.gz'


# In[894]:


# /usr/local/fsl/data/standard/


# In[895]:


# Testing

# res = MNI152_2mm_mask.run()
# res.outputs


# In[896]:


# afni.Bandpass.help()


# ## Band Pass Filtering
# Let's do a band pass filtering on the data using the code from https://neurostars.org/t/bandpass-filtering-different-outputs-from-fsl-and-nipype-custom-function/824/2

# In[897]:


### AFNI

bandpass = Node(afni.Bandpass(highpass=0.008, lowpass=0.08,
                         despike=False, no_detrend=True, notrans=True,
                         outputtype='NIFTI_GZ'),name='bandpass')

# bandpass = Node(afni.Bandpass(highpass=0.001, lowpass=0.01,
#                          despike=False, no_detrend=True, notrans=True,
#                          tr=2.0,outputtype='NIFTI_GZ'),name='bandpass')


# bandpass.inputs.mask = MNI152_2mm.outputs.mask_file


# In[898]:


# Testing bandpass on the func data in subject's space

# First comment out the bandpass.inputs.mask as it is in standard space.

# subject_id = layout.get_subjects()[0] # gives the first subject's ID
# func_file_path = [f.filename for f in layout.get(subject=subject_id, type='bold', extensions=['nii', 'nii.gz'])]
# bandpass.inputs.in_file = func_file_path[0]
# res = bandpass.run();


# In[899]:


# res.outputs.out_file


# In[900]:


# To view in fsl I need to save this file. You can change the the location as per your need.
# First run utility functions section. It contains the load_and_save function
# load_and_save(res.outputs.out_file,'/home/jovyan/work/preprocess/result/filtered_func.nii')


# In[901]:


# afni.Bandpass.help() # to see what all parameters are supported by Bandpass filter of afni


# ### Following is a Join Node that collects the preprocessed file paths and saves them in a file

# In[902]:


def save_file_list_function(in_brain, in_mask, in_motion_params, in_motion_outliers, in_joint_xformation_matrix, in_tr, in_atlas):
    # Imports
    import numpy as np
    import os
    from os.path import join as opj


    file_list = np.asarray(in_brain)
    print('######################## File List ######################: \n',file_list)

    np.save('brain_file_list',file_list)
    file_name = 'brain_file_list.npy'
    out_brain = opj(os.getcwd(),file_name) # path


    file_list2 = np.asarray(in_mask)
    print('######################## File List ######################: \n',file_list2)

    np.save('mask_file_list',file_list2)
    file_name2 = 'mask_file_list.npy'
    out_mask = opj(os.getcwd(),file_name2) # path


    file_list3 = np.asarray(in_motion_params)
    print('######################## File List ######################: \n',file_list3)

    np.save('motion_params_file_list',file_list3)
    file_name3 = 'motion_params_file_list.npy'
    out_motion_params = opj(os.getcwd(),file_name3) # path


    file_list4 = np.asarray(in_motion_outliers)
    print('######################## File List ######################: \n',file_list4)

    np.save('motion_outliers_file_list',file_list4)
    file_name4 = 'motion_outliers_file_list.npy'
    out_motion_outliers = opj(os.getcwd(),file_name4) # path


    file_list5 = np.asarray(in_joint_xformation_matrix)
    print('######################## File List ######################: \n',file_list5)

    np.save('joint_xformation_matrix_file_list',file_list5)
    file_name5 = 'joint_xformation_matrix_file_list.npy'
    out_joint_xformation_matrix = opj(os.getcwd(),file_name5) # path

    tr_list = np.asarray(in_tr)
    print('######################## TR List ######################: \n',tr_list)

    np.save('tr_list',tr_list)
    file_name6 = 'tr_list.npy'
    out_tr = opj(os.getcwd(),file_name6) # path


    file_list7 = np.asarray(in_atlas)
    print('######################## File List ######################: \n',file_list7)

    np.save('atlas_file_list',file_list7)
    file_name7 = 'atlas_file_list.npy'
    out_atlas = opj(os.getcwd(),file_name7) # path




    return out_brain, out_mask, out_motion_params, out_motion_outliers, out_joint_xformation_matrix, out_tr , out_atlas



# In[903]:


save_file_list = JoinNode(Function(function=save_file_list_function, input_names=['in_brain', 'in_mask', 'in_motion_params','in_motion_outliers','in_joint_xformation_matrix', 'in_tr', 'in_atlas'],
                 output_names=['out_brain','out_mask','out_motion_params','out_motion_outliers','out_joint_xformation_matrix','out_tr', 'out_atlas']),
                 joinsource="infosource",
                 joinfield=['in_brain', 'in_mask', 'in_motion_params','in_motion_outliers','in_joint_xformation_matrix','in_tr', 'in_atlas'],
                 name="save_file_list")


# In[904]:


# ------------------Change it in the program below -- all the names of parameters iin the workflow..


# ## AFNI's filter is working good:
# ### Next:
# - [x] Add the mask as parameter to the afni Node
# - [] Add the Node to the workflow
# - [x] Improve the data sink
# - [] Create Voxel pair FC map

# ### Motion outliers

# In[905]:


motionOutliers = Node(MotionOutliers(no_motion_correction=True, out_metric_plot = 'refrms_plot.png',
                                     out_metric_values='refrms_raw.txt'),name='motionOutliers')


# In[906]:


# (MotionOutliers(in_file = 'var.nii',no_motion_correction=False, out_metric_plot = 'refrms_plot',
#                                      out_metric_values='refrms_raw')).cmdline


# ### Lets see the number of regions in the atlas and display the atlas.

# In[907]:


# num_ROIs = int((np.max(atlas) - np.min(atlas) ))

# print('Min Index:', np.min(atlas),'Max Index', np.max(atlas))
# print('Total Number of Parcellations = ',num_ROIs)


# ## Workflow for atlas registration  from std to functional

# In[908]:


wf_atlas_resize_reg = Workflow(name=atlas_resize_reg_directory)

wf_atlas_resize_reg.connect([

            # Apply the inverse matrix to the 3mm Atlas to transform it to func space

            (maskfunc4mean, std2func_xform, [(('out_file','reference'))]),

            (resample_atlas, std2func_xform, [('out_file','in_file')] ),

            # Now, applying the inverse matrix

            (inv_mat, std2func_xform, [('out_file','in_matrix_file')]), # output: Atlas in func space

            (std2func_xform, save_file_list, [('out_file','in_atlas')]),

            # ---------------------------Save the required files --------------------------------------------

            (save_file_list, dataSink, [('out_motion_params','motion_params_paths.@out_motion_params')]),
            (save_file_list, dataSink, [('out_motion_outliers','motion_outliers_paths.@out_motion_outliers')]),
            (save_file_list, dataSink, [('out_brain','preprocessed_brain_paths.@out_brain')]),
            (save_file_list, dataSink, [('out_mask','preprocessed_mask_paths.@out_mask')]),

            (save_file_list, dataSink, [('out_joint_xformation_matrix',
                                         'joint_xformation_matrix_paths.@out_joint_xformation_matrix')]),

            (save_file_list, dataSink, [('out_tr','tr_paths.@out_tr')]),

            (save_file_list, dataSink, [('out_atlas','atlas_paths.@out_atlas')])



])


# In[909]:


wf_coreg_reg = Workflow(name=coreg_reg_directory)
# wf_coreg_reg.base_dir = base_directory
# Dir where all the outputs will be stored(inside coregistrationPipeline folder).


wf_coreg_reg.connect([


            # (BIDSDataGrabber,resample_anat,[('anat_file_path','in_file')]), # Resampled the anat file to 3mm
            #
            # (resample_anat,skullStrip,[('out_file','in_file')]),
            #
            # (skullStrip,func2anat_reg,[('out_file','reference')]), # Make the resampled file as reference in func2anat_reg
            #
            # # Sec 1. The above 3 steps registers the mean image to resampled anat image and
            # # calculates the xformation matrix .. I hope the xformation matrix will be saved
            #
            # (MNI152_2mm, resample_mni, [('standard_file','in_file')]),
            #
            # (resample_mni, anat2std_reg, [('out_file','reference')]),
            #
            # (skullStrip, anat2std_reg,  [('out_file','in_file')]),

            (BIDSDataGrabber,skullStrip,[('anat_file_path','in_file')]), # Resampled the anat file to 3mm

            (skullStrip,resample_anat,[('out_file','in_file')]),

            (resample_anat,func2anat_reg,[('out_file','reference')]), # Make the resampled file as reference in func2anat_reg

            # Sec 1. The above 3 steps registers the mean image to resampled anat image and
            # calculates the xformation matrix .. I hope the xformation matrix will be saved

            (MNI152_2mm, resample_mni, [('standard_file','in_file')]),

            (resample_mni, anat2std_reg, [('out_file','reference')]),

            (resample_anat, anat2std_reg,  [('out_file','in_file')]),


            # Calculates the Xformationmatrix from anat3mm to MNI 3mm

            # We can get those matrices by refering to func2anat_reg.outputs.out_matrix_file and similarly for anat2std_reg



            (func2anat_reg, concat_xform, [('out_matrix_file','in_file')]),

            (anat2std_reg, concat_xform, [('out_matrix_file','in_file2')]),

            (concat_xform, dataSink, [('out_file', 'tranformation_matrix_fun2std.@out_file')]),

            (concat_xform, save_file_list, [('out_file', 'in_joint_xformation_matrix')]),

            # Now inverse the func2std MAT to std2func
            (concat_xform, wf_atlas_resize_reg, [('out_file','inv_mat.in_file')])

# -------------------------------------------------



#             # Apply the transformation to 3mm Atlas to transform it to func space

#             (maskfunc4mean, std2func_xform, [(('standard_file','reference'))]),

#             (resample_atlas, std2func_xform, [('out_file','in_file')] ), # Applies the transform to the ...
#             # ... atlas

#             # Now, apply the inverse matrix to the atlas

#             (inv_mat, std2func_xform, [('out_file','in_matrix_file')]),


# ------------------------------------------



#             (save_file_list, dataSink, [('out_motion_params','motion_params_paths.@out_motion_params')]),
#             (save_file_list, dataSink, [('out_motion_outliers','motion_outliers_paths.@out_motion_outliers')]),
#             (save_file_list, dataSink, [('out_brain','preprocessed_brain_paths.@out_brain')]),
#             (save_file_list, dataSink, [('out_mask','preprocessed_mask_paths.@out_mask')]),

#             (save_file_list, dataSink, [('out_joint_xformation_matrix',
#                                          'joint_xformation_matrix_paths.@out_joint_xformation_matrix')]),

#             (save_file_list, dataSink, [('out_tr',
#                                          'tr_paths.@out_tr')])


])




# ## Co-Registration, Normalization and Bandpass Workflow
# 1. Co-registration means alligning the func to anat
# 2. Normalization means aligning func/anat to standard
# 3. Applied band pass filtering in range - highpass=0.008, lowpass=0.08

# In[910]:


wf_motion_correction_bet = Workflow(name=motion_correction_bet_directory)
# wf_motion_correction_bet.base_dir = base_directory

wf_motion_correction_bet.connect([

            (mcflirt,dataSink,[('par_file','motion_params.@par_file')]), # saves the motion parameters calculated before

            (mcflirt,save_file_list,[('par_file','in_motion_params')]),

#             (save_file_list, dataSink, [('out_motion_params','motion_params_paths.@out_motion_params')]),

            (mcflirt, meanfunc, [('out_file','in_file')]),

            (meanfunc, meanfuncmask, [('out_file','in_file')]),

            (slicetimer,motionOutliers,[('slice_time_corrected_file','in_file')]),

#             (mcflirt, motionOutliers, [('out_file','in_file')]),

            (meanfuncmask, motionOutliers, [('mask_file','mask')]),

            (motionOutliers, dataSink, [('out_file','motionOutliers.@out_file')]),

            (motionOutliers, dataSink, [('out_metric_plot','motionOutliers.@out_metric_plot')]),

            (motionOutliers, dataSink, [('out_metric_values','motionOutliers.@out_metric_values')]),

            (motionOutliers, save_file_list, [('out_file','in_motion_outliers')]),

#             (save_file_list, dataSink, [('out_motion_outliers','motion_outliers_paths.@out_motion_outliers')]),


            (mcflirt,applyMask , [('out_file','in_file')]), # 1

            (meanfuncmask, applyMask, [('mask_file','in_file2')]), # 2 output: 1&2,  BET on coregistered fmri scan

            (meanfunc, maskfunc4mean, [('out_file', 'in_file')]), # 3

            (meanfuncmask, maskfunc4mean, [('mask_file','in_file2')]), # 4 output: 3&4, BET on mean func scan


            (applyMask, save_file_list, [('out_file', 'in_brain')]),
            (applyMask, save_file_list, [('out_file2', 'in_mask')]),


#             (save_file_list, dataSink, [('out_brain','preprocessed_brain_paths.@out_brain')]),

#             (save_file_list, dataSink, [('out_mask','preprocessed_mask_paths.@out_mask')]),

            (maskfunc4mean, wf_coreg_reg, [('out_file','func2anat_reg.in_file')])



])


# ## Observation:
# Applying masking again on the Normalized func file greately reduced the size from ~600MB -> ~150MB. I think Normalizing might have generated some extra voxels in the region of 'no brain'. Masking again got rid of them. Hence, reduced size.

# ----------------------------------------------------------------------------------------------------------------------------
## Uncomment the below code to analyze only selected subjects to account for Slice time mismatch -----------------------------

## Adding a module to select only selected participants
#* KKI 
#* Leuven_1 
#* Leuven_2 
#* SBL 
#* Stanford (N.A)
#* Trinity 
#* UM_1 
#* UM_2 

import pandas as pd
#import numpy as np

df = pd.read_csv('/home1/varunk/data/ABIDE1/RawDataBIDs/composite_phenotypic_file.csv') # , index_col='SUB_ID'

df = df.sort_values(['SUB_ID'])
df
selected_participants = df.loc[(df['SITE_ID'] == 'KKI') | (df['SITE_ID'] == 'Leuven_1') | (df['SITE_ID'] == 'Leuven_2') \
                              | (df['SITE_ID'] == 'SBL') | (df['SITE_ID'] == 'Trinity') | (df['SITE_ID'] == 'UM_1') \
                              | (df['SITE_ID'] == 'UM_2')]

selected_participants = list(map(str, selected_participants.as_matrix(['SUB_ID']).squeeze()))

# ## Main Workflow

# In[911]:


subject_list = selected_participants#[0:2]
subject_list =  [str(item).zfill(7) for item in subject_list]

# ----------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------
# In[912]:


infosource = Node(IdentityInterface(fields=['subject_id']),
                  name="infosource")

infosource.iterables = [('subject_id',subject_list)]

# infosource.inputs.subject_id = subject_list[0]
# res = infosource.run()


# In[913]:


# res.outputs


# In[914]:


# Create the workflow
# Refer to Supplementary material's Section Two. for more on workspaces

wf = Workflow(name=parent_wf_directory)
# base_dir = opj(s,'result')
wf.base_dir = base_directory # Dir where all the outputs will be stored(inside BETFlow folder).

wf.connect([      (infosource, BIDSDataGrabber, [('subject_id','subject_id')]),
                  (BIDSDataGrabber, extract, [('func_file_path','in_file')]),
                  
                  (BIDSDataGrabber,getMetadata, [('func_file_path','in_file')]),
            
                  (getMetadata,slicetimer, [('tr','time_repetition')]),

                
                  (getMetadata,slicetimer, [('index_dir','index_dir')]),
            
                  (getMetadata,slicetimer, [('interleaved','interleaved')]),
            
                  (getMetadata,save_file_list, [('tr','in_tr')]),
                
                  (extract,slicetimer,[('roi_file','in_file')]),
                  (slicetimer,wf_motion_correction_bet,[('slice_time_corrected_file','mcflirt.in_file')])
           ])
# Run it in parallel
wf.run('MultiProc', plugin_args={'n_procs': 6})

# ### Summary:

# In[915]:


# Visualize the detailed graph
# from IPython.display import Image
wf.write_graph(graph2use='flat', format='png', simple_form=True)
# file_name = opj(base_directory,parent_wf_directory,'graph_detailed.dot.png')
# Image(filename=file_name)


# ### Summary [Incomplete]
# ```
# wf.connect([(infosource, BIDSDataGrabber, [('data_dir','data_dir'), ('subject_id', 'subject_id'),]),
#                   (BIDSDataGrabber, extract, [('func_file_path','in_file')]),
#                   (extract,slicetimer,[('roi_file','in_file')]),
#                   (slicetimer,mcflirt,[('slice_time_corrected_file','in_file')]),
#                   (mcflirt, skullStrip, [('mean_img', 'in_file')]),
#                   (mcflirt,applyMask,[('out_file','brain_file')]),
#                   (skullStrip, applyMask, [('mask_file', 'mask_file')]),
#                   ])
# ```
#
# In the above created workflow the `infosource` node iterates over the `subject_id`, it creates a Node and for each Subject ID it sends `data_dir` (path where the data resides) and the subject specific `subject_id` to `BIDSDataGrabber` Node.
#
# `BIDSDataGrabber` Node accepts the above 2 parameters, calls the function `get_nifti_filenames(subject_id,data_dir)`which returns the path of the anatomical and BOLD files of the subject with given subject_id and hence the Node produces output that I call `func_file_path` and `anat_file_path`. I have used only `func_file_path`right now.
#
# The file path denoted by '``func_file_path``' is then fed as input to `extract` that removes 4 initial brain volumes of the functional scan.
#
# Its output is called - `slice_time_corrected_file` which is fed to `mcflirt` node to correct the movion between volumes of an individual subject. This is called **Motion Correction**.
#
# In next step the mean_image from `mcflirt` is sent to `skullStrip` to get the mask. The role of `skullStrip` is just to obtain mask from the mean EPI image.
#
# The mask got above is then applied to the functional volume to get rif of skull.
#
#
#
# The final results are stored in the directory : `/home/jovyan/work/preprocess/result/BETFlow`. Every node has its own folder where its results are stored.
#
#
#

# # TODO
# * Make a single workflow
# * See how the file looks like after the transformation is applied
# * Change the FC code such that the FC maps is calculated instead of just matrices.
# * Transform these FC maps to Standard space
# * Do FDR correction (Check it with the MATLAB Fdr code)

# In[916]:


# os.chdir('../results/ABIDE1_Preprocess/motion_correction_bet/_subject_id_0050002/applyMask')

# in_file = 'sub-0050002_task-rest_run-1_bold_roi_st_mcf.nii_brain.nii.gz'

# in_matrix_file = 'sub-0050002_T1w_resample_brain_flirt_sub-0050002_task-rest_run-1_bold_roi_st_mcf_mean_bet_flirt.mat'



# In[917]:


# func2std_xform.inputs


# In[918]:

#
# import numpy as np
# X = np.load('../results_again/ABIDE1_Preprocess_Datasink/preprocessed_brain_paths/brain_file_list.npy')
#
#
# # In[919]:
#
#
# X
#
#
# # In[920]:
#
#
# X = np.load('../results_again/ABIDE1_Preprocess_Datasink/preprocessed_mask_paths/mask_file_list.npy')
# X
#
#
# # In[921]:
#
#
# get_ipython().system(' cat ../results_again/ABIDE1_Preprocess_Datasink/motion_params/sub-0050002/sub-0050002_motion_params.par')
#
#
# # In[922]:
#
#
# get_ipython().system(' cat ../results_again/ABIDE1_Preprocess_Datasink/tranformation_matrix_fun2std/sub-0050002/sub-0050002_task-rest_run-1_bold_roi_st_mcf_mean_bet_flirt_sub-0050002_T1w_resample_brain_flirt.mat')
#
#
# # In[923]:
#
#
# get_ipython().system(' cat ../results_again/ABIDE1_Preprocess_Datasink/tranformation_matrix_fun2std/sub-0050003/sub-0050003_task-rest_run-1_bold_roi_st_mcf_mean_bet_flirt_sub-0050003_T1w_resample_brain_flirt.mat')
#
#
# # $ \alpha_2 $ -- Just checking latex embedding
#
# # In[924]:
#
#
# #TR
# np.load('../results_again/ABIDE1_Preprocess_Datasink/tr_paths/tr_list.npy')
#
#
# # In[926]:
#
#
# np.load('../results_again/ABIDE1_Preprocess_Datasink/atlas_paths/atlas_file_list.npy')
#
#
# # In[927]:
#
#
# # Read the std atlas
#
# atlasObject = nib.load(atlasPath)
# atlas = atlasObject.get_data()
#
# num_ROIs = int((np.max(atlas) - np.min(atlas) ))
#
#
# # In[928]:
#
#
# num_ROIs
#
#
# # In[930]:
#
#
# # Read the func atlas
# atlaspath2 = '/home1/varunk/results_again/ABIDE1_Preprocess/motion_correction_bet/coreg_reg/atlas_resize_reg_directory/_subject_id_0050002/std2func_xform/fullbrain_atlas_thr0-2mm_resample_flirt.nii'
# atlasObject = nib.load(atlaspath2)
# atlas = atlasObject.get_data()
#
# num_ROIs = int((np.max(atlas) - np.min(atlas) ))
#
#
# # In[931]:
#
#
# num_ROIs
#
