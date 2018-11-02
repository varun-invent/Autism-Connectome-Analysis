import numpy as np

# Use the motion parameters to fid the bad brains
# Use the qc_csv to find the bad brains

#  Motion Based outliers


def read_par_file(motion_params_file):
    trans_x = []
    trans_y = []
    trans_z = []
    rot_x = []
    rot_y = []
    rot_z = []

    with open(motion_params_file) as f:
        for line in f:
            line = line.split(' ')
            # print(line)
            trans_x.append(float(line[6]))
            trans_y.append(float(line[8]))
            trans_z.append(float(line[10]))
            rot_x.append(float(line[0]))
            rot_y.append(float(line[2]))
            rot_z.append(float(line[4]))
    return trans_x, trans_y, trans_z, rot_x, rot_y, rot_z

def motion_outliers(motion_params_file, threshold):
    motion_params_paths = np.load(motion_params_npy)
    outliers_sub_ids = []
    for subject_param_path in motion_params_paths:
        sub_id = subject_param_path.split('/')[-1].split('_')[0].split('-')[1]
        trans_x, trans_y, trans_z, rot_x, rot_y, rot_z = read_par_file(subject_param_path)
        params = np.array([trans_x, trans_y, trans_z, rot_x, rot_y, rot_z]).T
        outlier_idx = np.where(np.abs(params) > thresh)
        num_outlier_entries = len(set(outlier_idx[0]))
        if num_outlier_entries >= 10:
            print('Subject %s with %s Outliers'%(sub_id, num_outlier_entries ))
            outliers_sub_ids.append(sub_id)
    return outliers_sub_ids

if __name__ == "__main__":
    motion_params_npy = '/mnt/project2/home/varunk/fMRI/results/resultsABIDE1_2/'+\
    'Preprocess_Datasink/motion_params_paths/motion_params_file_list.npy'

    thresh = 2.5
    outliers_sub_ids = motion_outliers(motion_params_npy, thresh)
    print(outliers_sub_ids)
        # print(num_outlier_entries)
