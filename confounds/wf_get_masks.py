from nipype.interfaces.fsl import (FLIRT, FAST, ConvertXFM, ImageMaths, MultiImageMaths, ApplyMask)
import nibabel as nib
import numpy as np
from scipy.ndimage.morphology import binary_erosion as erode
from nipype.pipeline import Node, Workflow
from nipype.interfaces.utility import IdentityInterface, Function
from nipype.interfaces.afni import Resample


def get_wf_tissue_masks(name='wf_tissue_masks'):
    '''
    This Function gives a workflow that resamples the T1 brains, extracts the
    tissue types thresholds at 0.5 and registers them to T2* space
    It then registers the tissue priors to the T2* space and then performs a
    bitwise AND between two maps.
    '''
    # csf_tissue_prior_path, gm_tissue_prior_path, wm_tissue_prior_path,
    # threshold = 0.5

    wf_tissue_masks = Workflow(name=name)

    inputspec = Node(IdentityInterface(
                fields=['resampled_anat_file_path', 'func2anat_mat_path','std2func_mat_path',
                'reference_func_file_path','brain_mask_eroded', 'threshold']),
                      name="inputspec")

    # FSL FAST node to segment the T1 brain
    fast = Node(FAST(out_basename = 'fast_'), name='fast')
    # probability_maps=True,segments=True,
    wf_tissue_masks.connect(inputspec, 'resampled_anat_file_path', fast, 'in_files')

    #  Invert the func2anat matrix to get anat2func
    inv_mat = Node(ConvertXFM(invert_xfm=True), name='inv_mat')
    wf_tissue_masks.connect(inputspec, 'func2anat_mat_path', inv_mat, 'in_file')

    # Transform the above segmented tissue masks to the functional space using the inverse matrix
    anat2func_xform_csf = Node(FLIRT(output_type='NIFTI',
                             apply_xfm=True, interp='sinc'), name='anat2func_xform_csf')

    wf_tissue_masks.connect(inputspec, 'reference_func_file_path', anat2func_xform_csf, 'reference')
    wf_tissue_masks.connect(inv_mat, 'out_file', anat2func_xform_csf, 'in_matrix_file')

    anat2func_xform_wm = Node(FLIRT(output_type='NIFTI',
                                apply_xfm=True, interp='sinc'), name='anat2func_xform_wm')
    wf_tissue_masks.connect(inputspec, 'reference_func_file_path', anat2func_xform_wm, 'reference')
    wf_tissue_masks.connect(inv_mat, 'out_file', anat2func_xform_wm, 'in_matrix_file')

    std2func_xform_eroded_brain = Node(FLIRT(output_type='NIFTI',
                                apply_xfm=True, interp='nearestneighbour'), name='std2func_xform_eroded_brain')
    wf_tissue_masks.connect(inputspec, 'reference_func_file_path', std2func_xform_eroded_brain, 'reference')
    wf_tissue_masks.connect(inputspec, 'std2func_mat_path', std2func_xform_eroded_brain, 'in_matrix_file')

    def select_item_from_array(arr, index=0):
        import numpy as np
        arr = np.array(arr)
        return arr[index]

    wf_tissue_masks.connect(fast, ('partial_volume_files', select_item_from_array, 0), anat2func_xform_csf, 'in_file')
    wf_tissue_masks.connect(fast, ('partial_volume_files', select_item_from_array, 2), anat2func_xform_wm, 'in_file')
    wf_tissue_masks.connect(inputspec, 'brain_mask_eroded', std2func_xform_eroded_brain, 'in_file')


    # Threshold

    def get_opstring(threshold):
        op = '-thr '+str(threshold)+' -bin'
        return op

    # print(inputspec.outputs)
    # ----- CSF ------

    threshold_csf = Node(interface=ImageMaths(suffix='_thresh'),
                       name='threshold_csf')
    # threshold_csf.inputs.op_string = '-thresh '+str(inputspec.outputs.threshold)+' -bin'
    wf_tissue_masks.connect(inputspec, ('threshold', get_opstring), threshold_csf, 'op_string' )
    wf_tissue_masks.connect(anat2func_xform_csf, 'out_file', threshold_csf, 'in_file')





    # ------- GM --------

    # threshold_gm = Node(interface=ImageMaths(op_string='-thresh',
    #                                             suffix='_thresh'),
    #                    name='threshold_gm')
    #
    #
    # wf_tissue_priors.connect(inputspec, ('threshold', get_opstring), threshold_gm, 'op_string' )
    # wf_tissue_priors.connect(fast, partial_volume_map[1], threshold_gm, 'in_file')
    #
    # -------- WM --------

    threshold_wm = Node(interface=ImageMaths(suffix='_thresh'),
                       name='threshold_wm')
    wf_tissue_masks.connect(inputspec, ('threshold', get_opstring), threshold_wm, 'op_string')
    wf_tissue_masks.connect(anat2func_xform_wm, 'out_file', threshold_wm, 'in_file')

    #  -------------------

    #
    # wf_tissue_masks.connect(threshold_csf, 'out_file', std2func_xform_csf, 'in_file')
    # wf_tissue_masks.connect(threshold_wm, 'out_file', std2func_xform_wm, 'in_file')


    # Masking the outer brain CSF

    csf_mask = Node(interface=ApplyMask(),
                       name='csf_mask')
    wf_tissue_masks.connect(threshold_csf, 'out_file', csf_mask, 'in_file')
    wf_tissue_masks.connect(std2func_xform_eroded_brain, 'out_file', csf_mask, 'mask_file')

    # Masking the outer brain WM that might be present due to poor BET

    wm_mask = Node(interface=ApplyMask(),
                       name='wm_mask')
    wf_tissue_masks.connect(threshold_wm, 'out_file', wm_mask, 'in_file')
    wf_tissue_masks.connect(std2func_xform_eroded_brain, 'out_file', wm_mask, 'mask_file')


    # wm_mask = Node(interface=ApplyMask(),
    #                    name='wm_mask')
    # wf_tissue_masks.connect(std2func_xform_wm, 'out_file', wm_mask, 'in_file')
    # wf_tissue_masks.connect(std2func_xform_wm_prior, 'out_file', wm_mask, 'mask_file')



    outputspec = Node(IdentityInterface(fields=['csf_mask', 'wm_mask']),
    name="outputspec")

    wf_tissue_masks.connect(csf_mask, 'out_file', outputspec, 'csf_mask')
    # wf_tissue_priors.connect(threshold_gm, 'out_file', outputspec, 'gm_tissue_prior_path')
    wf_tissue_masks.connect(wm_mask, 'out_file', outputspec, 'wm_mask')

    return wf_tissue_masks




if __name__ == "__main__":
    wf_tissue_masks = get_wf_tissue_masks(name='wf_tissue_masks3')

    wf_tissue_masks.inputs.inputspec.resampled_anat_file_path = \
    '/mnt/project1/home1/varunk/fMRI/testScripts/_subject_id_0050002/resample_anat/sub-0050002_T1w_brain_resample.nii'

    wf_tissue_masks.inputs.inputspec.reference_func_file_path = \
    '/mnt/project1/home1/varunk/fMRI/testScripts/func_subject_id_0050002/applyMask/sub-0050002_task-rest_run-1_bold_roi_st_mcf.nii_brain.nii.gz'

    wf_tissue_masks.inputs.inputspec.func2anat_mat_path = \
    '/mnt/project1/home1/varunk/fMRI/results/resultsABIDE1/preprocess/'+\
    'motion_correction_bet/coreg_reg/_subject_id_0050002/func2anat_reg/'+\
    'sub-0050002_task-rest_run-1_bold_roi_st_mcf_mean_bet_flirt.mat'

    wf_tissue_masks.inputs.inputspec.std2func_mat_path = \
    '/mnt/project1/home1/varunk/fMRI/results/resultsABIDE1/preprocess/'+\
    'motion_correction_bet/coreg_reg/atlas_resize_reg_directory/_subject_id_0050002/'+\
    'std2func_xform/fullbrain_atlas_thr0-2mm_resample_flirt.mat'


    # wf_tissue_masks.inputs.inputspec.tissue_prior_csf_path = '/mnt/project1/home1/varunk/fMRI/testScripts/results/wf_tissue_priors/threshold_csf/avg152T1_csf_resample_thresh.nii.gz'
    # wf_tissue_masks.inputs.inputspec.tissue_prior_wm_path = '/mnt/project1/home1/varunk/fMRI/testScripts/results/wf_tissue_priors/threshold_wm/avg152T1_white_resample_thresh.nii.gz'
    wf_tissue_masks.inputs.inputspec.brain_mask_eroded = '/mnt/project1/home1/varunk/fMRI/Autism-Connectome-Analysis/tissuepriors/brain_mask_2mm_eroded_18mm.nii.gz'

    wf_tissue_masks.inputs.inputspec.threshold = 0.2
    wf_tissue_masks.base_dir = 'results/'
    TEMP_DIR_FOR_STORAGE = 'crash_files/'
    wf_tissue_masks.config = {"execution": {"crashdump_dir": TEMP_DIR_FOR_STORAGE}}

    wf_tissue_masks.write_graph(graph2use='flat', format='png', simple_form=True)

    out = wf_tissue_masks.run()
