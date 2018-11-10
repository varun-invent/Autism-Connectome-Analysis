import nibabel as nib
import numpy as np
import os
from scipy.ndimage.morphology import binary_dilation as dilation
from scipy.ndimage.morphology import binary_erosion as erosion
import copy

# Take as input 3 tissue priors
# Merge them according to the highest probability
# Then seperate the 3 and convert them to masks and return them.

def merge(in_file_list, mask_file):
    priors_list = []
    affine = None
    header = None
    mask = nib.load(mask_file).get_data()
    for in_file in in_file_list:
        brain = nib.load(in_file)
        affine = brain.affine
        header = brain.header
        brain_data = brain.get_data().astype(np.float32)
        priors_list.append(brain_data)

    prior_list_np = np.array(priors_list)
    combined_prior = np.max(prior_list_np,axis=0)
    combined_prior_atlas = np.argmax(prior_list_np,axis=0)
    combined_prior_atlas[combined_prior != 0] = combined_prior_atlas[combined_prior != 0] + 1

    #  Apply mask on the combined_prior_atlas
    combined_prior_atlas = combined_prior_atlas * mask

    #  Seperate each of the priors

    out_files_list = []
    for i in range(1, len(in_file_list)+1):
        brain = copy.deepcopy(combined_prior_atlas)
        brain[combined_prior_atlas != i] = 0
        brain[combined_prior_atlas == i] = 1 # To create mask
        brain_file_name = 'prior_%s_thresh.nii.gz'%i
        out_file_brain = os.path.abspath(brain_file_name)

        brain_with_header = nib.Nifti1Image(brain, affine=affine,header = header)
        nib.save(brain_with_header,out_file_brain)
        out_files_list.append(out_file_brain)





    # brain_file_name = 'combined_tissue_priors_atlas.nii.gz'
    # out_file_brain = os.path.abspath(brain_file_name)
    #
    # brain_with_header = nib.Nifti1Image(combined_prior_atlas, affine=affine,header = header)
    # nib.save(brain_with_header,out_file_brain)

    # return out_file_brain
    return out_files_list

if __name__ ==  "__main__":
    base_path = '/mnt/project1/home1/varunk/fMRI/Autism-Connectome-Analysis/tissuepriors/'
    csf_prior = base_path + 'avg152T1_csf.nii.gz'
    wm_prior = base_path + 'avg152T1_white.nii.gz'
    gray_prior = base_path + 'avg152T1_gray.nii.gz'
    mask_file = base_path + 'MNI152_T1_2mm_brain_ero_dilM_16mm.nii.gz'

    res = merge([csf_prior,wm_prior ,gray_prior], mask_file)
