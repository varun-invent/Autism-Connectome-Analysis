import postProcessingFcVolMatchedModularDynamicPipeline as ppfc
# import preprocessingPipelineImprovedModular as prep
import preprocessingPipelineImprovedModularWOAnatDynamic as prep

import scoreCorrelation as sc

import hypothesisTest as ht
import fdrBrainResultsModular as fdres
import os
import json
from os.path import join as opj
import itertools
from bids.grabbids import BIDSLayout
import time
from pathlib import Path
import numpy as np
import shutil

# ------------- Paths ----------------------------------------------------------------------------------------

# New path to store results:
# results_path = '/mnt/scratch/varunk/'

# results_path = '/mnt/project1/home1/varunk/results/'

SELECT_SUBJECTS = False # Number of subjects to select from the whole randomly
number_of_selected_subjects = 2

# ----------------------------------------Don't Modify ------------------------------------------------------------


path_cwd = os.getcwd()
path_split_list = path_cwd.split('/')
s = path_split_list[0:-1] # for getting to the parent dir of pwd
s = opj('/',*s) # *s converts list to path, # very important to add '/' in the begining so it is read as directory later



json_path = 'scripts/json/paths.json'
with open(json_path, 'rt') as fp:
    task_info = json.load(fp)

base_directory = task_info["base_directory_for_results"]
motion_correction_bet_directory = task_info["motion_correction_bet_directory"]
parent_wf_directory = task_info["parent_wf_directory"]
functional_connectivity_directory = task_info["functional_connectivity_directory"]
# functional_connectivity_directory = 'temp_fc'
coreg_reg_directory = task_info["coreg_reg_directory"]
atlas_resize_reg_directory = task_info["atlas_resize_reg_directory"]

data_directory = task_info["data_directory"]

datasink_name = task_info["datasink_name"]
fc_datasink_name = task_info["fc_datasink_name"]
# fc_datasink_name = 'temp_dataSink'
atlasPath = task_info["atlas_path"] # Standard Brainnetome path


brain_path = opj(base_directory,datasink_name,'preprocessed_brain_paths/brain_file_list.npy')
mask_path = opj(base_directory,datasink_name,'preprocessed_mask_paths/mask_file_list.npy')
atlas_path = opj(base_directory,datasink_name,'atlas_paths/atlas_file_list.npy') # Brainetome atlas in functional space
tr_path = opj(base_directory,datasink_name,'tr_paths/tr_list.npy')
motion_params_path = opj(base_directory,datasink_name,'motion_params_paths/motion_params_file_list.npy')

func2std_mat_path = opj(base_directory, datasink_name,'joint_xformation_matrix_paths/joint_xformation_matrix_file_list.npy')

MNI3mm_path = opj(base_directory,parent_wf_directory,motion_correction_bet_directory,coreg_reg_directory,'resample_mni/MNI152_T1_2mm_brain_resample.nii')

hypothesis_test_dir = opj(base_directory, task_info["hypothesis_test_dir"])

fdr_results_dir = task_info["fdr_results_dir"]

score_corr_dir =  opj(base_directory,task_info["score_corr_dir"])

demographics_file_path = task_info['demographics_file_path']
phenotype_file_path = task_info['phenotype_file_path']
# categoryInfo = '/home1/varunk/data/NYU_Cocaine-BIDS/grouping.csv'
categoryInfo = None

# Tissue - CSF, WM and GM masks path
csf_path = opj(base_directory,datasink_name,'csf_mask_paths/csf_mask_file_list.npy')
wm_path = opj(base_directory,datasink_name,'wm_mask_paths/wm_mask_file_list.npy')

#  Binarized atlas mask path

binarized_atlas_mask_path = task_info["binarized_atlas_mask_path"]
OVERWRITE_POSTPROS_DIR = False


# --------------- Creating Log --------------------

# Get time

current_time = time.asctime( time.localtime(time.time()) )

if not os.path.exists(base_directory):
    os.makedirs(base_directory)

log_path = opj(base_directory,"log.txt")
log = open(log_path, 'a')
print("-------------Starting the analysis at %s--------------\n"%(current_time))
log.write("-------------Starting the analysis at %s--------------\n"%(current_time))
log.flush()



# ----------------------------------------------------------------------------------------------------------------------------------------------------------

vols = 120 # needed for volume matching. Tells how many volumes you need to keep.
number_of_skipped_volumes = 4 # TODO: Hardcoded inside preprocessing code: Make it dynamic in v2
#  True volumes will be vols - number_of_skipped_volumes


num_proc = 27


# number_of_subjects = -1
# number_of_subjects = 7 # Number of subjects you wish to work with

# ----------------------------- Getting Subjects -------------------------------
# ----------------------------------- BIDS -------------------------------------
layout = BIDSLayout(data_directory)


# if number_of_subjects == -1:
#     number_of_subjects = len(layout.get_subjects())

# ABIDE II Bugs - No Anat all, No func: 29622

bugs_abide2 = ['28093', '28093', '28681',  '28682', '28683',  '28687', '28711', '28712', '28713', '28741', '28745',  '28751', '28755', '28756', '28757', '28758',
'28759', '28761', '28762','28763', '28764','28765','28766','28767','28768','28769','28770','28771','28772','28773','28774','28775','28776','28777','28778','28779',
'28780','28781','28782','28783', '29622'
]


# ABIDE II Motion Outliers bugs with  >30% volumes >2.5mm or degree
bugs_abide2.extend(['28715', '28823', '29067', '29069', '29102', '29110', '29424', '29428', '29514', '29628', '29880', '29887', '29888', '29903'])


# Visual QC
bugs_abide2.extend(['30150', '30149', '29576', '29573', '29570', '29565', '29566', '29559', '29558',\
 '29557', '29554', '29547', '29546', '29545', '29541', '29538', '29056', '29055',\
  '29054', '29051', '29049', '29043', '29042', '29040', '29033', '29030', '29029',\
   '29028', '29027', '29026', '29025', '29021', '29020', '29019', '29016', '29015',\
    '28816', '28720', '28719', '28718', '28717', '28714', '28707', '28705', '28704',\
     '28701', '28698', '28695', '28693', '28691', '28689', '28686', '28685', '28684', '28676'])






# bugs_abide1 = ['51232','51233','51242','51243','51244','51245','51246','51247',
# '51270','51310','51276','50045', '50746', '50727', '50774',
# '0050313', '0051195']

bugs_abide1 = ['51232','51233','51242','51243','51244','51245','51246','51247',
'51270','51310', '50727']


'''
UCLA_1
'51232','51233' 51242','51243','51244',
'51245','51246', '51247', '51270', '51310': Anat not present, Update: Highres downloaded

Note: Highres looks like a T2 image. Not sure why the pixdim4 is equal to 5 in one of the subjects of UCLA_2

'51276' Functional file corrupt, Update: Downloaded new and copied

PITT
'50045' Looks fine to me use it again

Leuven_2
'50045', '50746' looks fine to me

'50727': Incomplete Brain Coverage

KKI
'50774': Functional brain file was corrupted while transfer. Downloaded a new one.

UM1:
'0050313': Movement issue

Stanford:
'0051195': Movements

So, The final bug in ABIDE I is 50727 of Leuven_2 due to partial brain coverage
Others - '51232','51233' 51242','51243','51244',
'51245','51246','51247','51270','51310' doesnot have anatomical file so were
discarded. So total of 11 subjects are discarded.
'''


bugs = bugs_abide2
# bugs = []
# subject_list = (layout.get_subjects())[0:number_of_subjects]
# subject_list = list(map(int, subject_list))
subject_list = layout.get_subjects()

# Ignore Bugs
subject_list = list(set(subject_list) - set(bugs))

# ----------------------------  Select 30 % of subjects ------------------
# import pdb; pdb.set_trace()
if SELECT_SUBJECTS == True:
    sub_list_filename = opj(base_directory,'sub_list.npy')


    if Path(sub_list_filename).exists():
        subject_list = np.load(sub_list_filename)
        print('Using the previously defined subject list')
    else:
        print('Creating a new subject list')
        index = np.arange(len(subject_list))
        np.random.shuffle(index)
        subject_list = subject_list[0:number_of_selected_subjects]
        np.save(sub_list_filename, subject_list)
# ------------------------------------------------------------------------
subject_list.sort()

print(subject_list)



# -----------------------------------File List----------------------------------
# group1FilesPath = ''
# group2FilesPath = ''
#
# group1FilesList = np.genfromtxt(group1FilesPath,dtype='unicode')
# group2FilesList = np.genfromtxt(group2FilesPath,dtype='unicode')
#
# fileList = group1FilesPath + group2FilesPath



# -----------------------------------------------------------------------------

paths = {
        'json_path' : json_path,
        'base_directory' : base_directory,
        'motion_correction_bet_directory' : motion_correction_bet_directory,
        'parent_wf_directory' : parent_wf_directory,
        'functional_connectivity_directory' : functional_connectivity_directory,
        'coreg_reg_directory' : coreg_reg_directory,
        'atlas_resize_reg_directory' : atlas_resize_reg_directory,
        'subject_list' : subject_list,
        'datasink_name' : datasink_name,
        'fc_datasink_name' : fc_datasink_name,
        'atlasPath' : atlasPath,
        'brain_path' : brain_path,
        'mask_path' : mask_path,
        'atlas_path' : atlas_path,
        'tr_path' : tr_path,
        'motion_params_path' : motion_params_path,
        'func2std_mat_path' : func2std_mat_path,
        'MNI3mm_path' : MNI3mm_path,
        'demographics_file_path' : demographics_file_path,
        'phenotype_file_path' : phenotype_file_path,
        'data_directory' : data_directory,
        'hypothesis_test_dir' : hypothesis_test_dir,
        'fdr_results_dir' : fdr_results_dir,
        'score_corr_dir' : score_corr_dir,
        'csf_path' : csf_path,
        'wm_path' : wm_path,
        'binarized_atlas_mask_path' : binarized_atlas_mask_path
}



PREPROC = 0
POSTPROC = 0
HYPOTEST = 1
FDRES = 1
SCORE_CORR = 0

# ABIDE II
run = 1
session = [1,2]

match = 1 # Age matching
applyFisher = True

# itr = (list(itertools.product([0, 1], repeat=3)))
# itr = [(1,1,1,1,1)]
# itr = [(1,0,0,0,0)]
# itr = [(1,1,1,1)]
itr = [(1,1,1,1)] # Post processing options

log.write("Operations:\n")
log.write("Preprocess = %s\n"%(PREPROC))
log.write("Postprocess = %s\n"%(POSTPROC))
log.write("Hypothesis Test = %s\n"%(HYPOTEST))
log.write("FDR correction and Vizualization = %s\n"%(FDRES))
log.write("Score-Connectivity Correlation  = %s\n"%(SCORE_CORR))
if SELECT_SUBJECTS == True:
    log.write("Working with only %s Subjects selected randomly\n "%(number_of_selected_subjects))
    # number_of_selected_subjects = 80
log.flush()

# ---------------------- Preprocess --------------------------------------------

ANAT = 1

# itr_preproc = [1,1,0,1]
itr_preproc = [1,1,0,1]
extract, slicetimer,motionOutliers, mcflirt= list(map(str, itr_preproc))
options_binary_string = extract+slicetimer+motionOutliers+mcflirt


# print('Calculating the BIDS data layout:')
# layout = BIDSLayout(data_directory) # TODO ( GIT Push)  Taken out as it takes lot of time to execute.


if PREPROC == 1:
    print('Preprocessing')

    log.write("Preprocessing Params\n")
    log.write("Remove begining slices  = %s\n"%(extract))
    log.write("Slice time correction  = %s\n"%(slicetimer))
    log.write("Calculate motionOutliers  = %s\n"%(motionOutliers)) # has a bug
    log.write("Do motion correction using McFLIRT  = %s\n"%(mcflirt))

    log.flush()

    DO_FAST = False
    try:
        prep.main(paths,options_binary_string, ANAT, DO_FAST, num_proc)
    except:
        print('Error Occured in Preprocessing')

    # try:
    # except shutil.SameFileError:
    #     pass


# Options for Calculating residuals
calc_residual_options = np.array(['csf', 'wm', 'motion', 'global']) # 'const'
residual_options_itr = list(itertools.product([True, False], repeat= 4))
# calc_residual_options_itr = [['const']]*len(residual_options_itr) # [['const'],['const'], ...]
# Above statement didnt work. It was creating copies of later inserted elements
calc_residual_options_itr = []
for i in range(len(residual_options_itr)):
    calc_residual_options_itr.append(['const'])
for i, mask in enumerate(residual_options_itr):
    calc_residual_options_itr[i].extend(calc_residual_options[list(mask)])

# ------------------------PostProcess------------------------------------------

'''
Brains that have > 30% volumes either translation or rotation > 2.5 mm or degree
Calculated using find_bad_brains.py script
'''
bugs_abide1.extend(['0050123', '0050279', '0050286',\
 '0050306', '0050489', '0051095', '0051213'])

'''
ABIDE 1 Bugs: (By manually looking at resampled and registered anatomical files)

50697 Neck
50694 Neck
50626 Neck
50625 Neck
50624 Neck
50618 Blurred and skull	Bad anat
50617 Neck
51324 Flipped brain	Due to improper skull strip
51296 Too blurred and spread and a lot of skull
51263 Flipped brain	Due to improper skull strip
50746 Flipped brain	Due to improper skull strip
'''

bugs_abide1.extend(['50697','50694','50626','50625','50624','50618','50617','51324','51296','51263','50746'])


bugs = bugs_abide2


if POSTPROC == 1 or HYPOTEST == 1 or FDRES == 1:
    print('PostProcessing')
    log.write("Postprocessing Params\n")

    print('Residual_options_list: ',calc_residual_options_itr)

    # overriding calc_residual_options_itr for testing
    # calc_residual_options_itr = [['const','csf', 'wm', 'motion', 'global']]
    # calc_residual_options_itr = [['const','csf', 'wm', 'global']]
    # calc_residual_options_itr = [['const']]

    # calc_residual_options_itr = [[]]

    for calc_residual_options in calc_residual_options_itr:
        for calc_residual, smoothing, band_pass_filtering, volCorrect in itr:
            log.write("calc_residual  = %s\n"%(calc_residual))
            log.write("calc_residual_options = %s \n"%(str(calc_residual_options)))
            # log.write("global_signal_regression  = %s\n"%(global_signal_regression))
            log.write("smoothing  = %s\n"%(smoothing))
            log.write("band_pass_filtering  = %s\n"%(band_pass_filtering))
            log.write("volCorrect  = %s\n"%(volCorrect))
            if volCorrect == 1:
                log.write("Vols for matching: %s\n"%(vols))
            log.write("Number_of_skipped_volumes: %s\n"%(number_of_skipped_volumes))
            log.write("Fisher Transform = %s\n"%(applyFisher))
            log.flush()

            #
            # comb = ''
            # for a in calc_residual_options:
            #     comb = comb + a
            #
            # combination = 'calc_residual' + str(int(calc_residual)) + \
            # 'smoothing' + str(int(smoothing)) +\
            # 'filt' + str(int(band_pass_filtering)) +\
            # 'calc_residual_options' + comb
            #
            # print("Combination: ",combination)
            # functional_connectivity_directory =  combination
            # print(calc_residual, smoothing,band_pass_filtering)
            save_npy = 0

            if POSTPROC == 1:
                ppfc.main(paths, vols, calc_residual, band_pass_filtering, smoothing, volCorrect, \
                number_of_skipped_volumes, num_proc, save_npy, calc_residual_options, OVERWRITE_POSTPROS_DIR,\
                run, session)


# ------------------- Hypothesis Test ------------------------------------------

            if HYPOTEST == 1:
                print('Hypothesis Test')
                log.write("Hypothesis Test\n")
                log.flush()

                # for calc_residual, smoothing,band_pass_filtering, volCorrect in itr:
                ht.main(paths, bugs,applyFisher,categoryInfo, match, calc_residual, band_pass_filtering, \
                    smoothing, num_proc, calc_residual_options, OVERWRITE_POSTPROS_DIR)


            # -------------------- FDR and results -----------------------------------------

            if FDRES == 1:
                print('FDR Correction and computing files for visualization of results')
                log.write("FDR Correction and computing files for visualization of results\n")
                log.flush()

                # for params in itr:
                fdres.main(paths, calc_residual, smoothing, band_pass_filtering,\
                 volCorrect, num_proc, calc_residual_options)



if SCORE_CORR == 1:
    print('Calculating Correlation-Score correaltions')

    # itr = (list(itertools.product([0, 1], repeat=3)))
    #
    # itr = [(1,1,0,1,1)]


    # bugs = []

    for calc_residual, smoothing,band_pass_filtering, volCorrect in itr:
        sc.main(paths, bugs,applyFisher,categoryInfo, match, calc_residual, band_pass_filtering, \
            smoothing, num_proc)
