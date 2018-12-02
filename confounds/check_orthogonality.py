import numpy as np
import nibabel as nib

#  Define Paths
confounds_path = '/mnt/project1/home1/varunk/fMRI/Autism-Connectome-Analysis/confounds/calc_residuals/'

mean_csf = confounds_path + 'mean_csf.txt'
mean_wm = confounds_path + 'mean_wm.txt'
mean_global = confounds_path + 'mean_global.txt'

residual_file = confounds_path + 'sub-0050009_task-rest_run-1_bold_roi_st_mcf.nii_brain.nii.gz_residual.nii.gz'

# Read The residual Brain using nib

residual_brain = nib.load(residual_file).get_data()

# Create a brain mask
global_mask = (residual_brain != 0).sum(-1) != 0
residual_brain_ndarray = residual_brain[global_mask] # Voxels x Time

# Read the csf, WM and Global mean using np

mean_csf_array = np.genfromtxt(mean_csf)
mean_wm_array = np.genfromtxt(mean_wm)
mean_global_array = np.genfromtxt(mean_global)


# Randomly sample 100 voxel time series by randomly selecting a list x from [0, brain.shape[0] - 1] and similarly y and z
sample_size = 20
voxel_list = np.random.choice(np.arange(residual_brain_ndarray.shape[0]),sample_size,replace=False)

# For these 100 voxels, calculate the dot product of the series and mean signal one by one
sum = 0
for voxel_idx in voxel_list:
            dot_prod = np.dot(residual_brain_ndarray[voxel_idx,:],mean_csf_array)
            sum = sum  + dot_prod
            print('Dot_prod: ',dot_prod)

print('Sum of Dot_prod ',sum)

#  CSF does not sum to be zero. Coz the orthogonalization was done using the all the confounds. THINK!
