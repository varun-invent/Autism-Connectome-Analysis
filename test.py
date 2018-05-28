# This file is for all kinds of Testing

# 1. Find how many subjects doesn't have anatomical files of run=1, session=1
# 2. Get the metadata for a given in_file. It is used to check if the metadata is being fetched up correctly or not.

from bids.grabbids import BIDSLayout

TEST = 1

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
