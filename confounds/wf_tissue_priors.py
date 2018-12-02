# from nipype.interfaces import fsl
from nipype.interfaces.fsl import (FLIRT, FAST, ConvertXFM, ImageMaths)
import nibabel as nib
import numpy as np
from scipy.ndimage.morphology import binary_erosion as erode
from nipype.pipeline import Node, Workflow
from nipype.interfaces.utility import IdentityInterface, Function
from nipype.interfaces.afni import Resample


def get_wf_tissue_priors(name='wf_tissue_priors3'):
    '''
    This Function gives a workflow that Resamples the tissue priors and then thresholds it at 0.5
    '''
    # csf_tissue_prior_path, gm_tissue_prior_path, wm_tissue_prior_path,
    # threshold = 0.5

    wf_tissue_priors = Workflow(name=name)

    inputspec = Node(IdentityInterface(fields=['csf_tissue_prior_path',  'wm_tissue_prior_path',
                                 'threshold','std2func_mat_path', 'reference_func_file_path']),
                      name="inputspec")
    '''
    # 'gm_tissue_prior_path',

    resample_tissue_prior_csf = Node(Resample(voxel_size=(3, 3, 3), resample_mode='Cu', # cubic interpolation
                             outputtype='NIFTI'),
                    name="resample_tissue_prior_csf")



    # resample_tissue_prior_gm = Node(Resample(voxel_size=(3, 3, 3), resample_mode='Cu', # cubic interpolation
    #                          outputtype='NIFTI'),
    #                 name="resample_tissue_prior_gm")



    resample_tissue_prior_wm = Node(Resample(voxel_size=(3, 3, 3), resample_mode='Cu', # cubic interpolation
                             outputtype='NIFTI'),
                    name="resample_tissue_prior_wm")


    wf_tissue_priors.connect(inputspec, 'csf_tissue_prior_path', resample_tissue_prior_csf, 'in_file' )
    # wf_tissue_priors.connect(inputspec, 'gm_tissue_prior_path', resample_tissue_prior_gm, 'in_file' )
    wf_tissue_priors.connect(inputspec, 'wm_tissue_prior_path', resample_tissue_prior_wm, 'in_file' )
    '''

    # #  Invert the func2anat matrix to get anat2func
    # inv_mat = Node(ConvertXFM(invert_xfm=True), name='inv_mat')
    # wf_tissue_priors.connect(inputspec, 'func2anat_mat_path', inv_mat, 'in_file')

    # Transform the  tissue priors to the functional space using the inverse matrix
    std2func_xform_csf_prior = Node(FLIRT(output_type='NIFTI',
                             apply_xfm=True, interp='sinc'), name='std2func_xform_csf_prior')

    wf_tissue_priors.connect(inputspec, 'reference_func_file_path', std2func_xform_csf_prior, 'reference')
    wf_tissue_priors.connect(inputspec, 'std2func_mat_path', std2func_xform_csf_prior, 'in_matrix_file')

    std2func_xform_wm_prior = Node(FLIRT(output_type='NIFTI',
                                apply_xfm=True, interp='sinc'), name='std2func_xform_wm_prior')
    wf_tissue_priors.connect(inputspec, 'reference_func_file_path', std2func_xform_wm_prior, 'reference')
    wf_tissue_priors.connect(inputspec, 'std2func_mat_path', std2func_xform_wm_prior, 'in_matrix_file')

    # Transformed the priors
    #  Get the input in_file(s) of the std2func_xform_csf and std2func_xform_wm from the old workspace
    wf_tissue_priors.connect(inputspec, 'csf_tissue_prior_path', std2func_xform_csf_prior, 'in_file')
    wf_tissue_priors.connect(inputspec, 'wm_tissue_prior_path', std2func_xform_wm_prior, 'in_file')





    # Threshold

    def get_opstring(threshold, tissue_type):
        if tissue_type == 'csf':
            max = 216  #  216 is the highest intensity of the resampled afni output for CSF
        elif tissue_type == 'wm':
            max = 253 #  253 is the highest intensity of the resampled afni output for WM

        threshold = int(threshold * max)
        op = '-thr '+str(threshold)+' -bin'
        return op

    # ----- CSF ------

    threshold_csf = Node(interface=ImageMaths(suffix='_thresh'),
                       name='threshold_csf')



    wf_tissue_priors.connect(inputspec, ('threshold', get_opstring, 'csf'), threshold_csf, 'op_string' )
    wf_tissue_priors.connect(std2func_xform_csf_prior, 'out_file', threshold_csf, 'in_file')

    # ------- GM --------

    # threshold_gm = Node(interface=ImageMaths(suffix='_thresh'),
    #                    name='threshold_gm')


    # wf_tissue_priors.connect(inputspec, ('threshold', get_opstring), threshold_gm, 'op_string' )
    # wf_tissue_priors.connect(resample_tissue_prior_gm, 'out_file', threshold_gm, 'in_file')

    # -------- WM --------

    threshold_wm = Node(interface=ImageMaths(suffix='_thresh'),
                       name='threshold_wm')

    wf_tissue_priors.connect(inputspec, ('threshold', get_opstring, 'wm'), threshold_wm, 'op_string' )
    wf_tissue_priors.connect(std2func_xform_wm_prior, 'out_file', threshold_wm, 'in_file')

    #  -------------------




    outputspec = Node(IdentityInterface(fields=['csf_tissue_prior_path', 'wm_tissue_prior_path', 'threshold']),
                      name="outputspec")

    # , 'gm_tissue_prior_path'
    wf_tissue_priors.connect(threshold_csf, 'out_file', outputspec, 'csf_tissue_prior_path')
    # wf_tissue_priors.connect(threshold_gm, 'out_file', outputspec, 'gm_tissue_prior_path')
    wf_tissue_priors.connect(threshold_wm, 'out_file', outputspec, 'wm_tissue_prior_path')

    return wf_tissue_priors

if __name__ == "__main__":
    tissue_priors = get_wf_tissue_priors()
    tissue_priors.inputs.inputspec.csf_tissue_prior_path = '/mnt/project1/home1/varunk/fMRI/Autism-Connectome-Analysis/tissuepriors/avg152T1_csf.nii.gz'
    # tissue_priors.inputs.inputspec.gm_tissue_prior_path = '/mnt/project1/home1/varunk/fMRI/Autism-Connectome-Analysis/tissuepriors/avg152T1_brain.nii.gz'
    tissue_priors.inputs.inputspec.wm_tissue_prior_path = '/mnt/project1/home1/varunk/fMRI/Autism-Connectome-Analysis/tissuepriors/avg152T1_white.nii.gz'
    tissue_priors.inputs.inputspec.threshold = 0.5
    # tissue_priors.inputs.inputspec.resampled_anat_file_path = \
    # '/mnt/project1/home1/varunk/fMRI/testScripts/_subject_id_0050002/resample_anat/sub-0050002_T1w_brain_resample.nii'

    tissue_priors.inputs.inputspec.reference_func_file_path = \
    '/mnt/project1/home1/varunk/fMRI/testScripts/func_subject_id_0050002/applyMask/sub-0050002_task-rest_run-1_bold_roi_st_mcf.nii_brain.nii.gz'

    tissue_priors.inputs.inputspec.std2func_mat_path = \
    '/mnt/project1/home1/varunk/fMRI/results/resultsABIDE1/preprocess/'+\
    'motion_correction_bet/coreg_reg/atlas_resize_reg_directory/_subject_id_0050002/'+\
    'std2func_xform/fullbrain_atlas_thr0-2mm_resample_flirt.mat'




    tissue_priors.base_dir = 'results/'
    TEMP_DIR_FOR_STORAGE = 'crash_files/'
    tissue_priors.config = {"execution": {"crashdump_dir": TEMP_DIR_FOR_STORAGE}}

    tissue_priors.write_graph(graph2use='flat', format='png', simple_form=True)

    out = tissue_priors.run()
