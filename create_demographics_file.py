import pandas as pd
import numpy as np
import json
import string
from bids.grabbids import BIDSLayout
import nibabel as nib

# --------------------------------------------------------------------- PART 1  GET Metadata-------------------------------------
'''
load BIDS data grabber
load the VOLUMES from metadata for each subject
Create a Sub_ID_VOlumes dict
select subjects that have volumes > threshold
'''


data_directory = '/mnt/project1/home1/varunk/data/ABIDE2RawDataBIDS'

layout = BIDSLayout(data_directory)
subjects = layout.get_subjects()

subid_vol_dict = {}
subject_list = []

session = [1,2]
run = 1
bugs_abide2 = ['28093', '28093', '28681',  '28682', '28683',  '28687', '28711', '28712', '28713', '28741', '28745',  '28751', '28755', '28756', '28757', '28758',
'28759', '28761', '28762','28763', '28764','28765','28766','28767','28768','28769','28770','28771','28772','28773','28774','28775','28776','28777','28778','28779',
'28780','28781','28782','28783', '29622'
]
for subject_id in subjects:
    func_file_path = [f.filename for f in layout.get(subject=subject_id, type='bold',session = session[0], run=run, extensions=['nii', 'nii.gz'])]
    if len(func_file_path) == 0:
        func_file_path = [f.filename for f in layout.get(subject=subject_id, type='bold',session = session[1], run=run, extensions=['nii', 'nii.gz'])]
        if len(func_file_path) == 0:
            if subject_id not in bugs_abide2:
                print('No Func file: %s'%subject_id)
            continue

    # print(func_file_path)
    metadata = layout.get_metadata(path=func_file_path[0])
    volumes  = metadata['NumberofMeasurements']
    try:
        volumes = int(volumes)
    except ValueError:
        # Mixed Volumes site
        brain_img = nib.load(func_file_path[0])
        volumes = brain_img.shape[-1]

    if volumes >= vols:
        subid_vol_dict[subject_id] = volumes
        subject_list.append(subject_id)


    print('Subject: %s Volumes: %s'%(subject_id, volumes))


import pdb; pdb.set_trace()



# df = pd.read_csv('/home1/varunk/data/ABIDE1/RawDataBIDs/composite_phenotypic_file.csv') # , index_col='SUB_ID'
phenotype_file_path = '/mnt/project1/home1/varunk/data/ABIDE2RawDataBIDS/ABIDEII_Composite_Phenotypic (copy).csv'
scan_params_file_path = '/mnt/project1/home1/varunk/data/ABIDE2RawDataBIDS/scan_params_file.txt'

df = pd.read_csv(phenotype_file_path, encoding = "ISO-8859-1")
df = df.sort_values(['SUB_ID'])

with open(scan_params_file_path, 'r') as f:
    scan_param_paths = f.read().split('\n')[0:-1]

print(scan_param_paths)

SITES = np.unique(df.as_matrix(['SITE_ID']).squeeze())
data_frame = pd.DataFrame({
'SITE_NAME': [] ,
'TR': [],
'VOLUMES': [],
'xdim_mm': [],
'ydim_mm': [],
'zdim_mm': [],
'xdim_voxels': [],
'ydim_voxels': [],
'zdim_voxels': [],
'NUM_AUT_DSM_V': [] ,
'NUM_AUT_MALE_DSM_V': [] ,
'NUM_AUT_FEMALE_DSM_V': [],
'NUM_AUT_AGE_lte12_DSM_V' : [],
'NUM_AUT_AGE_12_18_DSM_V' : [],
'NUM_AUT_AGE_18_24_DSM_V': [],
'NUM_AUT_AGE_24_34_DSM_V' :[],
'NUM_AUT_AGE_34_50_DSM_V' : [],
'NUM_AUT_AGE_gt50_DSM_V' : [],
'NUM_AUT_DSM_IV' : [],
'NUM_AUT_MALE_DSM_IV' : [],
'NUM_AUT_FEMALE_DSM_IV' : [],
'NUM_ASP_DSM_IV' : [],
'NUM_ASP_MALE_DSM_IV' : [],
'NUM_ASP_FEMALE_DSM_IV' : [],
'NUM_PDDNOS_DSM_IV' : [],
'NUM_PDDNOS_MALE_DSM_IV' : [],
'NUM_PDDNOS_FEMALE_DSM_IV' : [],
'NUM_ASP_PDDNOS_DSM_IV' : [],
'NUM_ASP_PDDNOS_MALE_DSM_IV' : [],
'NUM_ASP_PDDNOS_FEMALE_DSM_IV' : [],
'NUM_TD' : [],
'NUM_TD_MALE' : [],
'NUM_TD_FEMALE' : [],
'NUM_TD_AGE_lte12' : [],
'NUM_TD_AGE_12_18' : [],
'NUM_TD_AGE_18_24' : [],
'NUM_TD_AGE_24_34' : [],
'NUM_TD_AGE_34_50' : [],
'NUM_TD_AGE_gt50' : []

})

for SITE in SITES:
    NUM_AUT_DSM_V = df.loc[(df['DX_GROUP'] == 1) & (df['SITE_ID'] == SITE)].shape[0]
    NUM_AUT_MALE_DSM_V = df.loc[(df['DX_GROUP'] == 1) & (df['SEX'] == 1) & (df['SITE_ID'] == SITE)].shape[0]
    NUM_AUT_FEMALE_DSM_V = df.loc[(df['DX_GROUP'] == 1) & (df['SEX'] == 2) & (df['SITE_ID'] == SITE)].shape[0]

    NUM_AUT_AGE_lte12_DSM_V = df.loc[(df['DX_GROUP'] == 1) & (df['AGE_AT_SCAN'] <= 12) & (df['SITE_ID'] == SITE)].shape[0]
    NUM_AUT_AGE_12_18_DSM_V = df.loc[(df['DX_GROUP'] == 1) & (df['AGE_AT_SCAN'] > 12) & (df['AGE_AT_SCAN'] <= 18) & (df['SITE_ID'] == SITE)].shape[0]
    NUM_AUT_AGE_18_24_DSM_V = df.loc[(df['DX_GROUP'] == 1) & (df['AGE_AT_SCAN'] > 18) & (df['AGE_AT_SCAN'] <= 24) & (df['SITE_ID'] == SITE)].shape[0]
    NUM_AUT_AGE_24_34_DSM_V = df.loc[(df['DX_GROUP'] == 1) & (df['AGE_AT_SCAN'] > 24) & (df['AGE_AT_SCAN'] <= 34) & (df['SITE_ID'] == SITE)].shape[0]
    NUM_AUT_AGE_34_50_DSM_V = df.loc[(df['DX_GROUP'] == 1) & (df['AGE_AT_SCAN'] > 34) & (df['AGE_AT_SCAN'] <= 50) & (df['SITE_ID'] == SITE)].shape[0]
    NUM_AUT_AGE_gt50_DSM_V = df.loc[(df['DX_GROUP'] == 1) & (df['AGE_AT_SCAN'] > 50 ) & (df['SITE_ID'] == SITE)].shape[0]


    NUM_AUT_DSM_IV = df.loc[(df['DSM_IV_TR'] == 1) & (df['SITE_ID'] == SITE)].shape[0]
    NUM_AUT_MALE_DSM_IV = df.loc[(df['DSM_IV_TR'] == 1) & (df['SEX'] == 1) & (df['SITE_ID'] == SITE)].shape[0]
    NUM_AUT_FEMALE_DSM_IV = df.loc[(df['DSM_IV_TR'] == 1) & (df['SEX'] == 2) & (df['SITE_ID'] == SITE)].shape[0]


    NUM_ASP_DSM_IV = df.loc[(df['DSM_IV_TR'] == 2) & (df['SITE_ID'] == SITE)].shape[0]
    NUM_ASP_MALE_DSM_IV = df.loc[(df['DSM_IV_TR'] == 2) & (df['SEX'] == 1) & (df['SITE_ID'] == SITE)].shape[0]
    NUM_ASP_FEMALE_DSM_IV = df.loc[(df['DSM_IV_TR'] == 2) & (df['SEX'] == 2) & (df['SITE_ID'] == SITE)].shape[0]

    NUM_PDDNOS_DSM_IV = df.loc[(df['DSM_IV_TR'] == 3) & (df['SITE_ID'] == SITE)].shape[0]
    NUM_PDDNOS_MALE_DSM_IV = df.loc[(df['DSM_IV_TR'] == 3) & (df['SEX'] == 1) & (df['SITE_ID'] == SITE)].shape[0]
    NUM_PDDNOS_FEMALE_DSM_IV = df.loc[(df['DSM_IV_TR'] == 3) & (df['SEX'] == 2) & (df['SITE_ID'] == SITE)].shape[0]

    NUM_ASP_PDDNOS_DSM_IV = df.loc[(df['DSM_IV_TR'] == 4) & (df['SITE_ID'] == SITE)].shape[0]
    NUM_ASP_PDDNOS_MALE_DSM_IV = df.loc[(df['DSM_IV_TR'] == 4) & (df['SEX'] == 1) & (df['SITE_ID'] == SITE)].shape[0]
    NUM_ASP_PDDNOS_FEMALE_DSM_IV = df.loc[(df['DSM_IV_TR'] == 4) & (df['SEX'] == 2) & (df['SITE_ID'] == SITE)].shape[0]



    NUM_TD = df.loc[(df['DX_GROUP'] == 2) & (df['SITE_ID'] == SITE)].shape[0]
    NUM_TD_MALE = df.loc[(df['DX_GROUP'] == 2) & (df['SEX'] == 1) & (df['SITE_ID'] == SITE)].shape[0]
    NUM_TD_FEMALE = df.loc[(df['DX_GROUP'] == 2) & (df['SEX'] == 2) & (df['SITE_ID'] == SITE)].shape[0]

    NUM_TD_AGE_lte12 = df.loc[(df['DX_GROUP'] == 2) & (df['AGE_AT_SCAN'] <= 12) & (df['SITE_ID'] == SITE)].shape[0]
    NUM_TD_AGE_12_18 = df.loc[(df['DX_GROUP'] == 2) & (df['AGE_AT_SCAN'] > 12) & (df['AGE_AT_SCAN'] <= 18) & (df['SITE_ID'] == SITE)].shape[0]
    NUM_TD_AGE_18_24 = df.loc[(df['DX_GROUP'] == 2) & (df['AGE_AT_SCAN'] > 18) & (df['AGE_AT_SCAN'] <= 24) & (df['SITE_ID'] == SITE)].shape[0]
    NUM_TD_AGE_24_34 = df.loc[(df['DX_GROUP'] == 2) & (df['AGE_AT_SCAN'] > 24) & (df['AGE_AT_SCAN'] <= 34) & (df['SITE_ID'] == SITE)].shape[0]
    NUM_TD_AGE_34_50 = df.loc[(df['DX_GROUP'] == 2) & (df['AGE_AT_SCAN'] > 34) & (df['AGE_AT_SCAN'] <= 50) & (df['SITE_ID'] == SITE)].shape[0]
    NUM_TD_AGE_gt50 = df.loc[(df['DX_GROUP'] == 2) & (df['AGE_AT_SCAN'] > 50 ) & (df['SITE_ID'] == SITE)].shape[0]

    tr = 0
    volumes = 0
    xdim_mm = 0
    ydim_mm = 0
    zdim_mm = 0
    xdim_voxels = 0
    ydim_voxels = 0
    zdim_voxels = 0

    # Accessing scan details
    for json_path in scan_param_paths:
        extracted_site = json_path.split('/')[-2]

        if (SITE).lower() in (extracted_site).lower():
            with open(json_path, 'rt') as fp:
                print('Site matched with ',json_path)
                task_info = json.load(fp)

    # Accessing the contents:

                tr = task_info['RepetitionTime']
                volumes = task_info['NumberofMeasurements']
                try:
                    dim_mm = task_info['PixelSpacing'].split('x')
                    xdim_mm, ydim_mm = dim_mm[0], dim_mm[1]
                except IndexError:
                    dim_mm = task_info['PixelSpacing'].split('X')
                    xdim_mm, ydim_mm = dim_mm[0], dim_mm[1]

                zdim_mm = task_info['SpacingBetweenSlices']

                try:
                    xdim_voxels, ydim_voxels = task_info['AcquisitionMatrix'].split('x')
                except ValueError:
                    xdim_voxels, ydim_voxels = task_info['AcquisitionMatrix'].split('X')
                zdim_voxels = task_info['NumberOfSlices']


    _df = pd.DataFrame({
    'SITE_NAME': SITE ,
    'TR': tr ,
    'VOLUMES': volumes,
    'xdim_mm':xdim_mm,
    'ydim_mm':ydim_mm,
    'zdim_mm':zdim_mm,
    'xdim_voxels':xdim_voxels,
    'ydim_voxels':ydim_voxels,
    'zdim_voxels':zdim_voxels,

    'NUM_AUT_DSM_V': NUM_AUT_DSM_V ,
    'NUM_AUT_MALE_DSM_V': NUM_AUT_MALE_DSM_V ,
    'NUM_AUT_FEMALE_DSM_V': NUM_AUT_FEMALE_DSM_V,
    'NUM_AUT_AGE_lte12_DSM_V' : NUM_AUT_AGE_lte12_DSM_V,
    'NUM_AUT_AGE_12_18_DSM_V' : NUM_AUT_AGE_12_18_DSM_V,
    'NUM_AUT_AGE_18_24_DSM_V': NUM_AUT_AGE_18_24_DSM_V,
    'NUM_AUT_AGE_24_34_DSM_V' :NUM_AUT_AGE_24_34_DSM_V,
    'NUM_AUT_AGE_34_50_DSM_V' : NUM_AUT_AGE_34_50_DSM_V,
    'NUM_AUT_AGE_gt50_DSM_V' : NUM_AUT_AGE_gt50_DSM_V,
    'NUM_AUT_DSM_IV' : NUM_AUT_DSM_IV,
    'NUM_AUT_MALE_DSM_IV' : NUM_AUT_MALE_DSM_IV,
    'NUM_AUT_FEMALE_DSM_IV' : NUM_AUT_FEMALE_DSM_IV,
    'NUM_ASP_DSM_IV' : NUM_ASP_DSM_IV,
    'NUM_ASP_MALE_DSM_IV' : NUM_ASP_MALE_DSM_IV,
    'NUM_ASP_FEMALE_DSM_IV' : NUM_ASP_FEMALE_DSM_IV,
    'NUM_PDDNOS_DSM_IV' : NUM_PDDNOS_DSM_IV,
    'NUM_PDDNOS_MALE_DSM_IV' : NUM_PDDNOS_MALE_DSM_IV,
    'NUM_PDDNOS_FEMALE_DSM_IV' : NUM_PDDNOS_FEMALE_DSM_IV,
    'NUM_ASP_PDDNOS_DSM_IV' : NUM_ASP_PDDNOS_DSM_IV,
    'NUM_ASP_PDDNOS_MALE_DSM_IV' : NUM_ASP_PDDNOS_MALE_DSM_IV,
    'NUM_ASP_PDDNOS_FEMALE_DSM_IV' : NUM_ASP_PDDNOS_FEMALE_DSM_IV,
    'NUM_TD' : NUM_TD,
    'NUM_TD_MALE' : NUM_TD_MALE,
    'NUM_TD_FEMALE' : NUM_TD_FEMALE,
    'NUM_TD_AGE_lte12' : NUM_TD_AGE_lte12,
    'NUM_TD_AGE_12_18' : NUM_TD_AGE_12_18,
    'NUM_TD_AGE_18_24' : NUM_TD_AGE_18_24,
    'NUM_TD_AGE_24_34' : NUM_TD_AGE_24_34,
    'NUM_TD_AGE_34_50' : NUM_TD_AGE_34_50,
    'NUM_TD_AGE_gt50' : NUM_TD_AGE_gt50

    },index=[0],columns = [ 'SITE_NAME',
                            'TR',
                            'VOLUMES',
                            'xdim_mm',
                            'ydim_mm',
                            'zdim_mm',
                            'xdim_voxels',
                            'ydim_voxels',
                            'zdim_voxels',
                            'NUM_AUT_DSM_V',
                            'NUM_AUT_MALE_DSM_V',
                            'NUM_AUT_FEMALE_DSM_V',
                            'NUM_AUT_AGE_lte12_DSM_V',
                            'NUM_AUT_AGE_12_18_DSM_V',
                            'NUM_AUT_AGE_18_24_DSM_V',
                            'NUM_AUT_AGE_24_34_DSM_V',
                            'NUM_AUT_AGE_34_50_DSM_V',
                            'NUM_AUT_AGE_gt50_DSM_V',
                            'NUM_AUT_DSM_IV',
                            'NUM_AUT_MALE_DSM_IV',
                            'NUM_AUT_FEMALE_DSM_IV',
                            'NUM_ASP_DSM_IV',
                            'NUM_ASP_MALE_DSM_IV',
                            'NUM_ASP_FEMALE_DSM_IV',
                            'NUM_PDDNOS_DSM_IV',
                            'NUM_PDDNOS_MALE_DSM_IV',
                            'NUM_PDDNOS_FEMALE_DSM_IV',
                            'NUM_ASP_PDDNOS_DSM_IV',
                            'NUM_ASP_PDDNOS_MALE_DSM_IV',
                            'NUM_ASP_PDDNOS_FEMALE_DSM_IV',
                            'NUM_TD',
                            'NUM_TD_MALE',
                            'NUM_TD_FEMALE',
                            'NUM_TD_AGE_lte12',
                            'NUM_TD_AGE_12_18',
                            'NUM_TD_AGE_18_24',
                            'NUM_TD_AGE_24_34',
                            'NUM_TD_AGE_34_50',
                            'NUM_TD_AGE_gt50'])
    data_frame = data_frame.append(_df, ignore_index=True)[_df.columns.tolist()]

import pdb; pdb.set_trace()

data_frame.to_csv('demographics_ABIDE2.csv')
