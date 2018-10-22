import wf_tissue_priors as tp
import wf_get_masks as gm
from nipype.interfaces.fsl import (FLIRT, FAST, ConvertXFM, ImageMaths, MultiImageMaths)
import nibabel as nib
import numpy as np
from nipype.pipeline import Node, Workflow
from nipype.interfaces.utility import IdentityInterface, Function



#  Main Workflow that connects two workflows

def get_wf_main(name='wf_main'):

    wf_main = Workflow(name=name)

    inputspec = Node(IdentityInterface(fields=['resampled_anat_file_path',
                     'func2anat_mat_path', 'reference_func_file_path',
                      'csf_tissue_prior_path',  'wm_tissue_prior_path',
                      'threshold', 'std2func_mat_path', 'brain_mask_eroded']),
                      name="inputspec")
    outputspec = Node(IdentityInterface(fields=['csf_tissue_prior_path',  'wm_tissue_prior_path', 'qc_stats_dict']),
                      name="outputspec")

    tissue_priors = tp.get_wf_tissue_priors(name='wf_tissue_priors4')
    tissue_masks = gm.get_wf_tissue_masks(name='wf_tissue_masks4')


    def compute_qc_stats(anat_file_path, csf_mask, csf_prior, wm_mask, wm_prior):
        import numpy as np
        import nibabel as nib
        from collections import OrderedDict as od

        # print('$$$$$$$$$$$$Inside$$$$$$$$$QCFUNC')

        csf_prior_data = nib.load(csf_prior).get_data()

        wm_prior_data = nib.load(wm_prior).get_data()

        csf_mask_data = nib.load(csf_mask).get_data()

        wm_mask_data = nib.load(wm_mask).get_data()

        # A
        voxels_count_csf_prior = len((np.where(csf_prior_data == 1))[0])
        voxels_count_wm_prior = len((np.where(wm_prior_data == 1))[0])
        # B
        voxels_count_csf_mask = len((np.where(csf_mask_data == 1))[0])
        voxels_count_wm_mask = len((np.where(wm_mask_data == 1))[0])
        # A - B
        A_minus_B_csf =  len(np.where((csf_prior_data - csf_mask_data) == 1)[0])
        A_minus_B_wm = len(np.where((wm_prior_data - wm_mask_data) == 1)[0])
        # B - A
        B_minus_A_csf =  len(np.where((csf_prior_data - csf_mask_data) == -1)[0])
        B_minus_A_wm = len(np.where((wm_prior_data - wm_mask_data) == -1)[0])
        # A U B
        A_union_B_csf =  len(np.where((csf_prior_data + csf_mask_data) != 0)[0])
        A_union_B_wm = len(np.where((wm_prior_data + wm_mask_data) != 0)[0])
        # A I B
        A_intersection_B_csf =  len(np.where((csf_prior_data * csf_mask_data) == 1)[0])
        A_intersection_B_wm = len(np.where((wm_prior_data * wm_mask_data) == 1)[0])

        print('voxels_count_csf_prior ',voxels_count_csf_prior)
        print('voxels_count_wm_prior ',voxels_count_wm_prior)
        print('voxels_count_csf_mask ',voxels_count_csf_mask)
        print('voxels_count_wm_mask ',voxels_count_wm_mask)
        print('A_minus_B_csf ',A_minus_B_csf)
        print('A_minus_B_wm ',A_minus_B_wm)
        print('B_minus_A_csf ',B_minus_A_csf)
        print('B_minus_A_wm ',B_minus_A_wm)
        print('A_union_B_csf ',A_union_B_csf)
        print('A_union_B_wm ',A_union_B_wm)
        print('A_intersection_B_csf ',A_intersection_B_csf)
        print('A_intersection_B_wm ',A_intersection_B_wm)
        quality_csf =  A_intersection_B_csf/A_union_B_csf
        quality_wm =  A_intersection_B_wm/A_union_B_wm
        print('quality_csf ',quality_csf)
        print('quality_wm ',quality_wm)

        print('Anat File path ',anat_file_path)
        sub_id = anat_file_path.split('/')[-1].split('_')[0].split('-')[1]
        print('Sub ID ',sub_id)
        dict = od()
        dict['voxels_count_csf_prior'] = voxels_count_csf_prior
        dict['voxels_count_wm_prior'] = voxels_count_wm_prior
        dict['voxels_count_csf_mask'] = voxels_count_csf_mask
        dict['voxels_count_wm_mask'] = voxels_count_wm_mask
        dict['A_minus_B_csf'] = A_minus_B_csf
        dict['A_minus_B_wm'] = A_minus_B_wm
        dict['B_minus_A_csf'] = B_minus_A_csf
        dict['B_minus_A_wm'] = B_minus_A_wm
        dict['A_union_B_csf'] = A_union_B_csf
        dict['A_union_B_wm'] = A_union_B_wm
        dict['A_intersection_B_csf'] =  A_intersection_B_csf
        dict['A_intersection_B_wm'] = A_intersection_B_wm
        dict['quality_csf'] = quality_csf
        dict['quality_wm'] = quality_wm
        dict['sub_id'] = sub_id

        return dict


    qc_stats = Node(Function(function=compute_qc_stats, input_names=['anat_file_path',
     'csf_mask','csf_prior','wm_mask', 'wm_prior'],
    output_names=['dict']), name='qc_stats')


    wf_main.connect(inputspec, 'csf_tissue_prior_path',tissue_priors,'inputspec.csf_tissue_prior_path')
    wf_main.connect(inputspec, 'wm_tissue_prior_path',tissue_priors,'inputspec.wm_tissue_prior_path')
    wf_main.connect(inputspec, 'threshold',tissue_priors,'inputspec.threshold')
    wf_main.connect(inputspec, 'reference_func_file_path',tissue_priors,'inputspec.reference_func_file_path')
    wf_main.connect(inputspec, 'std2func_mat_path',tissue_priors,'inputspec.std2func_mat_path')

    wf_main.connect(tissue_priors, 'outputspec.csf_tissue_prior_path',outputspec,'csf_tissue_prior_path')
    wf_main.connect(tissue_priors, 'outputspec.wm_tissue_prior_path',outputspec,'wm_tissue_prior_path')


    wf_main.connect(inputspec, 'resampled_anat_file_path',tissue_masks,'inputspec.resampled_anat_file_path')
    wf_main.connect(inputspec, 'reference_func_file_path',tissue_masks,'inputspec.reference_func_file_path')
    wf_main.connect(inputspec, 'func2anat_mat_path',tissue_masks,'inputspec.func2anat_mat_path')
    wf_main.connect(inputspec, 'std2func_mat_path',tissue_masks,'inputspec.std2func_mat_path')
    wf_main.connect(inputspec, 'brain_mask_eroded',tissue_masks,'inputspec.brain_mask_eroded')
    wf_main.connect(inputspec, 'threshold',tissue_masks,'inputspec.threshold')

    wf_main.connect(tissue_masks, 'outputspec.csf_mask', outputspec,'csf_mask')
    wf_main.connect(tissue_masks, 'outputspec.wm_mask', outputspec,'wm_mask')

    wf_main.connect(tissue_priors, 'outputspec.csf_tissue_prior_path',qc_stats,'csf_prior')
    wf_main.connect(tissue_priors, 'outputspec.wm_tissue_prior_path',qc_stats,'wm_prior')
    wf_main.connect(tissue_masks, 'outputspec.csf_mask', qc_stats,'csf_mask')
    wf_main.connect(tissue_masks, 'outputspec.wm_mask', qc_stats,'wm_mask')
    wf_main.connect(inputspec, 'resampled_anat_file_path', qc_stats, 'anat_file_path')

    wf_main.connect(qc_stats, 'dict', outputspec, 'qc_stats_dict')

    return wf_main

#
# def func_qc_stats_to_csv(dict):
    # pass
#
# qc_stats_to_csv = JoinNode(Function(function=func_qc_stats_to_csv, input_names=['dict'],
#                output_names=['csv_path']),
#                joinsource="infosource",
#                joinfield=['in_brain'],
#                name="save_file_list_in_brain")


# --------------------------------------------------------------------------------------------------------------

wf_main = get_wf_main(name='wf_main')

wf_main.inputs.inputspec.resampled_anat_file_path = \
'/mnt/project1/home1/varunk/fMRI/testScripts/_subject_id_0050002/resample_anat/sub-0050002_T1w_brain_resample.nii'

wf_main.inputs.inputspec.reference_func_file_path = \
'/mnt/project1/home1/varunk/fMRI/testScripts/func_subject_id_0050002/applyMask/sub-0050002_task-rest_run-1_bold_roi_st_mcf.nii_brain.nii.gz'

wf_main.inputs.inputspec.func2anat_mat_path = \
'/mnt/project1/home1/varunk/fMRI/results/resultsABIDE1/preprocess/'+\
'motion_correction_bet/coreg_reg/_subject_id_0050002/func2anat_reg/'+\
'sub-0050002_task-rest_run-1_bold_roi_st_mcf_mean_bet_flirt.mat'

wf_main.inputs.inputspec.std2func_mat_path = \
'/mnt/project1/home1/varunk/fMRI/results/resultsABIDE1/preprocess/'+\
'motion_correction_bet/coreg_reg/atlas_resize_reg_directory/_subject_id_0050002/'+\
'std2func_xform/fullbrain_atlas_thr0-2mm_resample_flirt.mat'


# wf_tissue_masks.inputs.inputspec.tissue_prior_csf_path = '/mnt/project1/home1/varunk/fMRI/testScripts/results/wf_tissue_priors/threshold_csf/avg152T1_csf_resample_thresh.nii.gz'
# wf_tissue_masks.inputs.inputspec.tissue_prior_wm_path = '/mnt/project1/home1/varunk/fMRI/testScripts/results/wf_tissue_priors/threshold_wm/avg152T1_white_resample_thresh.nii.gz'
wf_main.inputs.inputspec.brain_mask_eroded = \
'/mnt/project1/home1/varunk/fMRI/Autism-Connectome-Analysis/tissuepriors/brain_mask_2mm_eroded_18mm.nii.gz'

wf_main.inputs.inputspec.threshold = 0.5

wf_main.inputs.inputspec.csf_tissue_prior_path =\
'/mnt/project1/home1/varunk/fMRI/Autism-Connectome-Analysis/tissuepriors/avg152T1_csf.nii.gz'
wf_main.inputs.inputspec.wm_tissue_prior_path =\
'/mnt/project1/home1/varunk/fMRI/Autism-Connectome-Analysis/tissuepriors/avg152T1_white.nii.gz'


wf_main.base_dir = 'results/'
TEMP_DIR_FOR_STORAGE = 'crash_files/'
wf_main.config = {"execution": {"crashdump_dir": TEMP_DIR_FOR_STORAGE}}

wf_main.write_graph(graph2use='flat', format='png', simple_form=True)

out = wf_main.run()
