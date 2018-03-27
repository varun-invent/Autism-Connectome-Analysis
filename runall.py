



import postProcessingFcVolMatchedModularDynamicPipeline as ppfc
# import preprocessingPipelineImprovedModular as prep
import preprocessingPipelineImprovedModularWOAnatDynamic as prep

import hypothesisTest as ht
import fdrBrainResultsModular as fdres
import os
import json
from os.path import join as opj
import itertools
from bids.grabbids import BIDSLayout

# ------------- Paths ----------------------------------------------------------------------------------------

path_cwd = os.getcwd()
path_split_list = path_cwd.split('/')
s = path_split_list[0:-1] # for getting to the parent dir of pwd
s = opj('/',*s) # *s converts list to path, # very important to add '/' in the begining so it is read as directory later



json_path = 'scripts/json/paths.json'
with open(json_path, 'rt') as fp:
    task_info = json.load(fp)



base_directory = opj(s,task_info["base_directory_for_results"])
motion_correction_bet_directory = task_info["motion_correction_bet_directory"]
parent_wf_directory = task_info["parent_wf_directory"]
functional_connectivity_directory = task_info["functional_connectivity_directory"]
# functional_connectivity_directory = 'temp_fc'
coreg_reg_directory = task_info["coreg_reg_directory"]
atlas_resize_reg_directory = task_info["atlas_resize_reg_directory"]
data_directory = opj(s,task_info["data_directory"])
datasink_name = task_info["datasink_name"]
fc_datasink_name = task_info["fc_datasink_name"]
# fc_datasink_name = 'temp_dataSink'
atlasPath = opj(s,task_info["atlas_path"]) # Standard Brainnetome path


brain_path = opj(base_directory,datasink_name,'preprocessed_brain_paths/brain_file_list.npy')
mask_path = opj(base_directory,datasink_name,'preprocessed_mask_paths/mask_file_list.npy')
atlas_path = opj(base_directory,datasink_name,'atlas_paths/atlas_file_list.npy') # Brainetome atlas in functional space
tr_path = opj(base_directory,datasink_name,'tr_paths/tr_list.npy')
motion_params_path = opj(base_directory,datasink_name,'motion_params_paths/motion_params_file_list.npy')

func2std_mat_path = opj(base_directory, datasink_name,'joint_xformation_matrix_paths/joint_xformation_matrix_file_list.npy')

MNI3mm_path = opj(base_directory,parent_wf_directory,motion_correction_bet_directory,coreg_reg_directory,'resample_mni/MNI152_T1_2mm_brain_resample.nii')

hypothesis_test_dir = opj(base_directory, task_info["hypothesis_test_dir"])

fdr_results_dir = task_info["fdr_results_dir"]

demographics_file_path = '/home1/varunk/Autism-Connectome-Analysis-brain_connectivity/notebooks/demographics.csv'
phenotype_file_path = '/home1/varunk/data/ABIDE1/RawDataBIDs/composite_phenotypic_file.csv'
# categoryInfo = '/home1/varunk/data/NYU_Cocaine-BIDS/grouping.csv'
categoryInfo = None

# ----------------------------------------------------------------------------------------------------------------------------------------------------------

vols = 120 # needed for volume matching. Tells how many volumes you need to keep.
number_of_skipped_volumes = 4
#  True volumes will be vols - number_of_skipped_volumes


num_proc = 4


number_of_subjects = -1
# number_of_subjects = 2 # Number of subjects you wish to work with



# ----------------------------- Getting Subjects -------------------------------
# ----------------------------------- BIDS -------------------------------------
layout = BIDSLayout(data_directory)


if number_of_subjects == -1:
    number_of_subjects = len(layout.get_subjects())


subject_list = (layout.get_subjects())[0:number_of_subjects]
# subject_list = list(map(int, subject_list))

# -----------------------------------File List----------------------------------
# group1FilesPath = ''
# group2FilesPath = ''
#
# group1FilesList = np.genfromtxt(group1FilesPath,dtype='unicode')
# group2FilesList = np.genfromtxt(group2FilesPath,dtype='unicode')
#
# fileList = group1FilesPath + group2FilesPath



# -----------------------------------------------------------------------------

paths = [json_path,
base_directory,
motion_correction_bet_directory,
parent_wf_directory,
functional_connectivity_directory,
coreg_reg_directory,
atlas_resize_reg_directory,
subject_list,
datasink_name,
fc_datasink_name,
atlasPath,
brain_path,
mask_path,
atlas_path,
tr_path,
motion_params_path,
func2std_mat_path,
MNI3mm_path,
demographics_file_path,
phenotype_file_path,
data_directory,
hypothesis_test_dir,
fdr_results_dir]

PREPROC = 1
POSTPROC = 1
HYPOTEST = 1
FDRES = 1

match = 1 # Age matching
applyFisher = True

# itr = (list(itertools.product([0, 1], repeat=3)))
itr = [(1,0,1,1,1)]
# itr = [(1,1,1,1,1)]
# ---------------------- Preprocess --------------------------------------------

ANAT = 1

itr_preproc = [1,1,0,1]
# itr_preproc = [0,0,0]
extract, slicetimer,motionOutliers, mcflirt= list(map(str, itr_preproc))
options_binary_string = extract+slicetimer+motionOutliers+mcflirt
if PREPROC == 1:
    print('Preprocessing')
    prep.main(paths,options_binary_string, ANAT, num_proc)


# ------------------------PostProcess------------------------------------------
if POSTPROC == 1:
    print('PostProcessing')



    for motion_param_regression, global_signal_regression, smoothing, band_pass_filtering, volCorrect in itr:
        combination = 'motionRegress' + str(int(motion_param_regression)) + \
         'global' + str(int(global_signal_regression)) + 'smoothing' + str(int(smoothing)) +\
         'filt' + str(int(band_pass_filtering))

        print("Combination: ",combination)
        functional_connectivity_directory =  combination
        print(motion_param_regression,  global_signal_regression, smoothing,band_pass_filtering)
        ppfc.main(paths, vols, motion_param_regression, global_signal_regression, band_pass_filtering, smoothing, volCorrect, \
        number_of_skipped_volumes, num_proc)


# ------------------- Hypothesis Test ------------------------------------------

if HYPOTEST == 1:
    print('Hypothesis Test')

    # itr = (list(itertools.product([0, 1], repeat=3)))
    #
    # itr = [(1,1,0,1,1)]

    bugs = ['51232','51233','51242','51243','51244','51245','51246','51247','51270','51310','50045', '51276', '50746', '50727', '51276']

    # bugs = []

    for motion_param_regression, global_signal_regression, smoothing,band_pass_filtering, volCorrect in itr:
        ht.main(paths, bugs,applyFisher,categoryInfo, match, motion_param_regression, global_signal_regression, band_pass_filtering, \
            smoothing, num_proc)

# -------------------- FDR and results -----------------------------------------

if FDRES == 1:
    print('FDR Correction and computing files for visualization of results')
    # itr = (list(itertools.product([0, 1], repeat=3)))

    # itr = [(1,1,0,1,1)]
    for params in itr:
        fdres.main(paths, params, num_proc = 7)
