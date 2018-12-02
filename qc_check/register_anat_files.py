import nibabel as nib
import numpy as np
from nipype.pipeline import Node, MapNode, Workflow
from nipype.interfaces.utility import IdentityInterface, Function
from nipype.interfaces.fsl import FLIRT
import os

# Node for applying xformation matrix to functional data

def wf_transform_anat(in_file_list, in_matrix_file_list, reference):
    func2std_xform = MapNode(FLIRT(output_type='NIFTI',
                             apply_xfm=True), name="func2std_xform",
                             iterfield = ['in_file', 'in_matrix_file', 'reference'])

    inputspec = Node(IdentityInterface(fields=['in_file_list', 'in_matrix_file_list', 'reference']),
                      name="inputspec")

    inputspec.inputs.in_file_list = in_file_list
    inputspec.inputs.in_matrix_file_list = in_matrix_file_list
    inputspec.inputs.reference = reference

    wf_transform_anat = Workflow(name="wf_transform_anat")
    wf_transform_anat.connect(inputspec,'in_file_list' , func2std_xform,'in_file')
    wf_transform_anat.connect(inputspec,'in_matrix_file_list' , func2std_xform,'in_matrix_file')
    wf_transform_anat.connect(inputspec,'reference' , func2std_xform,'reference')

    return wf_transform_anat


if __name__ == "__main__":
    in_file_list_path = '/mnt/project2/home/varunk/fMRI/results/'+\
    'resultsABIDE1_2/Preprocess_Datasink/resample_anat_filenames.txt'
    in_matrix_file_list_path = '/mnt/project2/home/varunk/fMRI/results/'+\
    'resultsABIDE1_2/Preprocess_Datasink/anat2std_mat_files.txt'
    reference_path = '/mnt/project2/home/varunk/fMRI/results/resultsABIDE1_2/'+\
    'preprocess/motion_correction_bet/coreg_reg/resample_mni/MNI152_T1_2mm_brain_resample.nii'

    in_file_list = [i.decode('UTF-8') for i in np.genfromtxt(in_file_list_path, dtype=None)]
    in_matrix_file_list = [i.decode('UTF-8') for i in np.genfromtxt(in_matrix_file_list_path, dtype=None)]
    reference = [reference_path]*len(in_matrix_file_list)

    print(in_file_list,'\n',in_matrix_file_list, '\n', reference)

    wf = wf_transform_anat(in_file_list, in_matrix_file_list, reference)
    wf.base_dir = '/mnt/project2/home/varunk/fMRI/results/'
    wf.run()
