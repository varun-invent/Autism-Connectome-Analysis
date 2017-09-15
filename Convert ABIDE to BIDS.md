# Conversion of ABIDE's NYU site data to BIDS format (Ubuntu)

The following tutorial is specific to data from ABIDE's NYU site but if understood nicely, can be extended to all INDI datasets such as ABIDE I, ABIDE II, CORR, ADHD200, ACPI and Rockland.


I downloaded only the NYU's site data from ABIDE and stored it in folder named `ABIDE`. To see its directory structure, in the terminal, I wrote :

```
$tree ABIDE
```
Result :

```
ABIDE
└── NYU
    ├── 0050952
    │   └── session_1
    │       ├── anat_1
    │       │   └── mprage.nii.gz
    │       └── rest_1
    │           └── rest.nii.gz
    ├── 0050953
    │   └── session_1
    │       ├── anat_1
    │       │   └── mprage.nii.gz
    │       └── rest_1
    │           └── rest.nii.gz
    ├── 0050954
    │   └── session_1
    │       ├── anat_1
    │       │   └── mprage.nii.gz
    │       └── rest_1
    │           └── rest.nii.gz
    ├── 0050955
    │   └── session_1
    │       ├── anat_1
    │       │   └── mprage.nii.gz
    │       └── rest_1
    │           └── rest.nii.gz
    ...
```
As you can see, it is similar to BIDS but not BIDS.

### Steps to convert it to BIDS format.

1. The github link of repository [indi_bidsification](https://github.com/FCP-INDI/indi_bidsification) contains scripts related to the process of converting INDI datasets into BIDS format. You will need the following files:
  - `indi_bidsification/ABIDE/participants_tsvs/nyu/participants.tsv`  
  This file contains details of all the participants such as Subject ID, Age, Gender, Handedness etc.
  - `indi_bidsification/ABIDE/scan_jsons/NYU/T1w.json`  
  This file contains detail related to anatomical MRI scan such as Acquision Time, Echo Time, Pixel Spacing etc.
  - `indi_bidsification/ABIDE/scan_jsons/NYU/task-rest_bold.json`  
  This file contains detail related to functional MRI scan such as Acquision Duration, Slice Thickness, Slice Timing etc.
  - `bidsify.sh`
  This file contains commands to actually convert ABIDE to BIDS format. While executing, I changed the name of input and output folder.  
    - In order to execute the commands in this file you need to download another python script called `pathmatcher.py` and it's dependencies form the github repository [neuro_experiments_tools](https://github.com/lrq3000/neuro_experiments_tools/tree/master/pathmatcher).
    - You can check the documentation on how to use `pathmatcher.py` [here](https://github.com/lrq3000/neuro_experiments_tools#regular-expression-path-matcher).  
- Download the whole folder in your active directory that contains the `ABIDE` folder.

- Create a new folder `ABIDE-BIDS`  

  ```
  $mkdir ABIDE-BIDS
  ```
- Change directory to the folder that contains `pathmatcher.py` and it's dependencies.
```
$cd pathmatcher/
```
- Execute the following commands:
```
$python pathmatcher.py -ri "(.+)/(\d{7})/session_(\d)/anat_(\d)/mprage.nii.gz" -ro "\1/sub-\2/anat/sub-\2_T1w.nii.gz" -i ABIDE -o ABIDE-BIDS -c
$python pathmatcher.py -ri "(.+)/(\d{7})/session_(\d)/rest_(\d)/rest.nii.gz" -ro "\1/sub-\2/func/sub-\2_task-rest_run-\4_bold.nii.gz" -i ../ABIDE -o ../ABIDE-BIDS -c
```
As you can see, the above commands are exactly same as the ones given in `bidsify.py` with just the folder names changed.  
The first command converts the folders containing the anatomical MRI files to BIDS and second command works on functional MRI files and folders.   

- Copy manually the following files:
  - `indi_bidsification/ABIDE/participants_tsvs/nyu/participants.tsv`  
  - `indi_bidsification/ABIDE/scan_jsons/NYU/T1w.json`  
  - `indi_bidsification/ABIDE/scan_jsons/NYU/task-rest_bold.json`  

which you downloaded from `indi_bidsification` to the NYU folder of ABIDE-BIDS created by the above mentioned `pathmatcher` commands.

- To check if BIDS was successfully created, execute the `tree` command in terminal.
  ```
  $tree ABIDE-BIDS
  ```
  This was the result:
  ```
  ABIDE-BIDS
└── NYU
    ├── participants.tsv
    ├── sub-0050952
    │   ├── anat
    │   │   └── sub-0050952_T1w.nii.gz
    │   └── func
    │       └── sub-0050952_task-rest_run-1_bold.nii.gz
    ├── sub-0050953
    │   ├── anat
    │   │   └── sub-0050953_T1w.nii.gz
    │   └── func
    │       └── sub-0050953_task-rest_run-1_bold.nii.gz
    ├── sub-0050954
    │   ├── anat
    │   │   └── sub-0050954_T1w.nii.gz
    │   └── func
    │       └── sub-0050954_task-rest_run-1_bold.nii.gz
    ...
    ├── T1w.json
    └── task-rest_bold.json

  ```

  So, now you know how to convert the NYU's site data of ABIDE to BIDS format. I you followed carefully, you can extend this procedure for converting other INDI datasets to BIDS. 

  Happy Learning!
