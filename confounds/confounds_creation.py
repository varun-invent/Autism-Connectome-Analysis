from nipype import Node
from nipype.interfaces.utility import Function

def calc_residuals(in_file, motion_file, csf_mask_path, gm_mask_path, wm_mask_path, global_signal_flag= None):
    import nibabel as nb
    import numpy as np
    import os
    from os.path import join as opj
    nii = nb.load(in_file)
    data = nii.get_data().astype(np.float32)
    global_mask = (data != 0).sum(-1) != 0

    # Calculate regressors
    regressor_map = {'constant' : np.ones((data.shape[3],1))}


    # Check and define regressors which are provided from files
    if motion_file is not None:
        motion = np.genfromtxt(motion_file)
        if motion.shape[0] != data.shape[3]:
            raise ValueError('Motion parameters {0} do not match data '
                             'timepoints {1}'.format(motion.shape[0],
                                                     data.shape[3]))
        if motion.size == 0:
            raise ValueError('Motion signal file {0} is '
                             'empty'.format(motion_file))
        regressor_map['motion'] = motion

    #  CSF
    if csf_mask_path is not None:
        csf_mask = nb.load(csf_mask_path).get_data().astype(np.float32)
        idx = np.where(csf_mask != 0)
        csf_signals = data[idx[0],idx[1],idx[2],:]
        # This should give me a matrix of size (T x CSF_Voxels)
        # Now take a mean of the above created matrix to get  a column vector i.e axis = 1

        mean_csf = np.mean(csf_signals.T, axis=1)

        if mean_csf.shape[0] != data.shape[3]:
            raise ValueError('CSF parameters {0} do not match data '
                             'timepoints {1}'.format(mean_csf.shape[0],
                                                     data.shape[3]))

        regressor_map['csf'] = mean_csf

    #  GM
    # if gm_mask_path is not None:
    #     gm_mask = nb.load(gm_mask_path).get_data().astype(np.float32)
    #     idx = np.where(gm_mask != 0)
    #     gm_signals = data[idx[0],idx[1],idx[2],:]
    #     # This should give me a matrix of size (T x gm_Voxels)
    #     # Now take a mean of the above created matrix to get  a column vector i.e axis = 1
    #
    #     mean_gm = np.mean(gm_signals.T, axis=1)
    #
    #     if mean_gm.shape[0] != data.shape[3]:
    #         raise ValueError('GM parameters {0} do not match data '
    #                          'timepoints {1}'.format(mean_gm.shape[0],
    #                                                  data.shape[3]))



    #  WM
    if wm_mask_path is not None:
        wm_mask = nb.load(wm_mask_path).get_data().astype(np.float32)
        idx = np.where(wm_mask != 0)
        wm_signals = data[idx[0],idx[1],idx[2],:]
        # This should give me a matrix of size (T x wm_Voxels)
        # Now take a mean of the above created matrix to get  a column vector i.e axis = 1

        mean_wm = np.mean(wm_signals.T, axis=1)

        if mean_wm.shape[0] != data.shape[3]:
            raise ValueError('wm parameters {0} do not match data '
                             'timepoints {1}'.format(mean_wm.shape[0],
                                                     data.shape[3]))

        regressor_map['wm'] = mean_wm

    if global_signal_flag is not None:
        idx = np.where(global_mask != 0)
        global_signal = data[idx[0],idx[1],idx[2],:]
        # This should give me a matrix of size (T x wm_Voxels)
        # Now take a mean of the above created matrix to get  a column vector i.e axis = 1

        mean_global = np.mean(global_signal.T, axis=1)
        regressor_map['global'] = mean_global

        print('Global_signal Shape: ',global_signal.shape)






    X = np.zeros((data.shape[3], 1))

    for rname, rval in regressor_map.items():
        X = np.hstack((X, rval.reshape(rval.shape[0],-1)))

    X = X[:,1:]

    if np.isnan(X).any() or np.isnan(X).any():
        raise ValueError('Regressor file contains NaN')

    Y = data[global_mask].T

    try:
        B = np.linalg.inv(X.T.dot(X)).dot(X.T).dot(Y)
    except np.linalg.LinAlgError as e:
        if "Singular matrix" in e:
            raise Exception("Error details: {0}\n\nSingular matrix error: "
                            "The nuisance regression configuration you "
                            "selected may have been too stringent, and the "
                            "regression could not be completed. Ensure your "
                            "parameters are not too "
                            "extreme.\n\n".format(e))
        else:
            raise Exception("Error details: {0}\n\nSomething went wrong with "
                            "nuisance regression.\n\n".format(e))

    Y_res = Y - X.dot(B)
    print('Shape of Y is: ', Y.shape)
    print('Shape of X is: ', X.shape)
    print('Shape of B is: ', B.shape)

    data[global_mask] = Y_res.T

    img = nb.Nifti1Image(data, header=nii.get_header(),
                         affine=nii.get_affine())

    subject_name = in_file
    # .split('/')[-1].split('.')[0]
    filename = subject_name + '_residual.nii.gz'
    out_file = os.path.join(os.getcwd(),filename )
    img.to_filename(out_file) # alt to nib.save

    return out_file

if __name__ == "__main__":
    calc_residuals = Node(Function(function=calc_residuals, input_names=['in_file', 'motion_file', 'csf_mask_path', 'gm_mask_path', 'wm_mask_path','global_signal_flag'],
                                    output_names=['out_file']), name='calc_residuals')


    calc_residuals.inputs.in_file = '/mnt/project1/home1/varunk/fMRI/results/temp_resultsABIDE1/preprocess/motion_correction_bet/_subject_id_0050009/applyMask/sub-0050009_task-rest_run-1_bold_roi_st_mcf.nii_brain.nii.gz'
    calc_residuals.inputs.motion_file = '/mnt/project1/home1/varunk/fMRI/results/temp_resultsABIDE1/preprocess/_subject_id_0050009/mcflirt/sub-0050009_task-rest_run-1_bold_roi_st_mcf.nii.par'
    calc_residuals.inputs.csf_mask_path = '/mnt/project1/home1/varunk/fMRI/results/temp_resultsABIDE1/preprocess/motion_correction_bet/coreg_reg/_subject_id_0050009/fast/fast__seg_0_flirt.nii'
    calc_residuals.inputs.gm_mask_path = '/mnt/project1/home1/varunk/fMRI/results/temp_resultsABIDE1/preprocess/motion_correction_bet/coreg_reg/_subject_id_0050009/fast/fast__seg_1_flirt.nii'
    calc_residuals.inputs.wm_mask_path = '/mnt/project1/home1/varunk/fMRI/results/temp_resultsABIDE1/preprocess/motion_correction_bet/coreg_reg/_subject_id_0050009/fast/fast__seg_2_flirt.nii'
    calc_residuals.inputs.global_signal_flag = True
    calc_residuals.run()
