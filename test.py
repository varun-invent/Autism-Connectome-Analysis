# This file is for all kinds of Testing

# 1. Find how many subjects doesn't have anatomical files of run=1, session=1
# 2. Get the metadata for a given in_file. It is used to check if the metadata is being fetched up correctly or not.
# 3. Checks if the corresponding brain npy files are same or not.
# 4. Checks if the corresponding brain nii files are same or not.

from bids.grabbids import BIDSLayout
import numpy as np
import nibabel as nib
from tqdm import tqdm

TEST = 4

################################################################################################################
                                             # Test 1
################################################################################################################

if TEST == 1:
    # from tqdm.tqdm import tqdm

    data_dir = "/mnt/project1/home1/varunk/data/ABIDE2RawDataBIDS"
    layout = BIDSLayout(data_dir)

    subject_list = layout.get_subjects()
    print('Subject List:')
    # print(subject_list)

    def get_nifti_filenames(subject_id):
        run = 1
        session = 1

        anat_file_path = [f.filename for f in layout.get(subject=subject_id, type='T1w', session = session, run=run, extensions=['nii', 'nii.gz'])]
        func_file_path = [f.filename for f in layout.get(subject=subject_id, type='bold',session = session, run=run, extensions=['nii', 'nii.gz'])]


        if len(anat_file_path) == 0:
            print('Anatomical file not present for subject ID %s' % subject_id )
        if len(func_file_path)  == 0:
            print('Functional file not present for subject ID %s' % subject_id )



    # Run
    for sub in subject_list:
        get_nifti_filenames(sub)


################################################################################################################
                                             # Test 2
################################################################################################################

if TEST == 2:

    # in_file = '/mnt/project1/home1/varunk/data/ABIDE2RawDataBIDS/ABIDEII-KKI_1/sub-29274/ses-1/func/sub-29274_ses-1_task-rest_acq-rc8chan_run-1_bold.nii.gz'

    in_file = '/mnt/project1/home1/varunk/data/ABIDE2RawDataBIDS/ABIDEII-KKI_1/sub-29474/ses-1/func/sub-29474_ses-1_task-rest_acq-rc32chan_run-1_bold.nii.gz'

    data_directory = '/mnt/project1/home1/varunk/data/ABIDE2RawDataBIDS'
    layout = BIDSLayout(data_directory)
    metadata = layout.get_metadata(path=in_file)

    print('Metadata:')
    print(metadata)


################################################################################################################
                                             # Test 3
################################################################################################################

if TEST == 3:
    in_file1 = np.load('../resultsABIDE_verifyCode/fc_datasink/motionRegress1global0smoothing1filt1/fc_map_npy_file_list.npy')
    in_file2 = np.load('../resultsABIDE_verifyCode2/fc_datasink/motionRegress1global0smoothing1filt1/fc_map_npy_file_list.npy')

    # in_file1 = np.load('../resultsABIDE_verifyCode/preprocess/_subject_id_0050774/mcflirt/sub-0050774_task-rest_run-1_bold_roi_st_mcf.nii')
    # in_file2 = np.load('../resultsABIDE_verifyCode2/preprocess/_subject_id_0050774/mcflirt/sub-0050774_task-rest_run-1_bold_roi_st_mcf.nii')

    # in_file1 = '../resultsABIDE_verifyCode/hypothesis_test/motionRegress1global0smoothing1filt1/Tvals.npy'
    # in_file2 = '../resultsABIDE_verifyCode2/hypothesis_test/motionRegress1global0smoothing1filt1/Tvals.npy'

    # volume = None

    # import pdb;pdb.set_trace()
    volume = 20
    slice = 30

    if isinstance(in_file1, np.ndarray):
        for file1, file2 in zip(in_file1,in_file2):
            # print('File 1: %s and File 2: %s'%(file1,file2))
            # brain1 = nib.load(file1).get_data()
            # brain2 = nib.load(file2).get_data()

            brain1 = np.load(file1,mmap_mode = 'r+')
            brain2 = np.load(file2,mmap_mode = 'r+')

            # import pdb; pdb.set_trace()
            brain1_slice = np.diag(brain1[:,:,slice,volume])
            brain2_slice = np.diag(brain2[:,:,slice,volume])

            import pdb;pdb.set_trace()
            if (brain1_slice == brain2_slice).all() :
                print('All slices are same')
            else:
                print('Brain Slices different - %s %s'%(file1,file2))

    else:
        brain1 = np.load(in_file1, mmap_mode = 'r+')
        brain2 = np.load(in_file2, mmap_mode = 'r+')

        brain1_slice = np.diag(brain1[:,:,slice,volume])
        brain2_slice = np.diag(brain2[:,:,slice,volume])

        if (brain1_slice == brain2_slice).all() :
            print('All slices are same')
        else:
            print('Brain Slices different - %s %s'%(in_file1,in_file2))



################################################################################################################
                                             # Test 4
################################################################################################################

if TEST == 4:
    # in_file1 = np.load('../resultsABIDE_verifyCode/fc_datasink/motionRegress1global0smoothing1filt1/fc_map_npy_file_list.npy')
    # in_file2 = np.load('../resultsABIDE_verifyCode2/fc_datasink/motionRegress1global0smoothing1filt1/fc_map_npy_file_list.npy')

    in_file1 = '../resultsABIDE_verifyCode/preprocess/_subject_id_0051277/mcflirt/sub-0051277_task-rest_run-1_bold_roi_st_mcf.nii'
    in_file2 = '../resultsABIDE_verifyCode2/preprocess/_subject_id_0051277/mcflirt/sub-0051277_task-rest_run-1_bold_roi_st_mcf.nii'

    # in_file1 = '../resultsABIDE_verifyCode/hypothesis_test/motionRegress1global0smoothing1filt1/Tvals.npy'
    # in_file2 = '../resultsABIDE_verifyCode2/hypothesis_test/motionRegress1global0smoothing1filt1/Tvals.npy'

    # volume = None

    # import pdb;pdb.set_trace()
    volume = 20
    slice = 30

    if isinstance(in_file1, np.ndarray):
        for file1, file2 in zip(in_file1,in_file2):
            # print('File 1: %s and File 2: %s'%(file1,file2))
            # brain1 = nib.load(file1).get_data()
            # brain2 = nib.load(file2).get_data()

            brain1 = nib.load(file1).get_data()
            brain2 = nib.load(file2).get_data()

            # import pdb; pdb.set_trace()
            brain1_slice = np.diag(brain1[:,:,slice,volume])
            brain2_slice = np.diag(brain2[:,:,slice,volume])

            import pdb;pdb.set_trace()
            if (brain1_slice == brain2_slice).all() :
                print('All slices are same')
            else:
                print('Brain Slices different - %s %s'%(file1,file2))

    else:
        brain1 = nib.load(in_file1).get_data()
        brain2 = nib.load(in_file2).get_data()

        brain1_slice = np.diag(brain1[:,:,slice,volume])
        brain2_slice = np.diag(brain2[:,:,slice,volume])

        if (brain1_slice == brain2_slice).all() :
            print('All slices are same')
        else:
            print('Brain Slices different - %s %s'%(in_file1,in_file2))

























'''

bugs_abide2 = ['28093', '28093', '28681',  '28682', '28683',  '28687', '28711', '28712', '28713', '28741', '28745',  '28751', '28755', '28756', '28757', '28758',
'28759', '28761', '28762','28763', '28764','28765','28766','28767','28768','28769','28770','28771','28772','28773','28774','28775','28776','28777','28778','28779',
'28780','28781','28782','28783','29622'
]


Anatomical file not present for subject ID 28093
Functional file not present for subject ID 28093
Functional file not present for subject ID 28681
Anatomical file not present for subject ID 28682
Functional file not present for subject ID 28683
Functional file not present for subject ID 28687
Functional file not present for subject ID 28711
Functional file not present for subject ID 28712
Functional file not present for subject ID 28713
Functional file not present for subject ID 28741
Functional file not present for subject ID 28745
Functional file not present for subject ID 28751
Functional file not present for subject ID 28755
Functional file not present for subject ID 28756
Functional file not present for subject ID 28757
Functional file not present for subject ID 28758
Functional file not present for subject ID 28759
Functional file not present for subject ID 28761
Functional file not present for subject ID 28762
Functional file not present for subject ID 28763
Functional file not present for subject ID 28764
Functional file not present for subject ID 28765
Functional file not present for subject ID 28766
Functional file not present for subject ID 28767
Functional file not present for subject ID 28768
Functional file not present for subject ID 28769
Functional file not present for subject ID 28770
Functional file not present for subject ID 28771
Functional file not present for subject ID 28772
Functional file not present for subject ID 28773
Functional file not present for subject ID 28774
Functional file not present for subject ID 28775
Functional file not present for subject ID 28776
Functional file not present for subject ID 28777
Functional file not present for subject ID 28778
Functional file not present for subject ID 28779
Functional file not present for subject ID 28780
Functional file not present for subject ID 28781
Functional file not present for subject ID 28782
Functional file not present for subject ID 28783
Functional file not present for subject ID 29622
'''
