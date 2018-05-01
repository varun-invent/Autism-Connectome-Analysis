# Issues with ABIDE-II data.

## Issue 1

I downloaded the ABIDE-II data from Amazon bucket using cyberduck on Mac. There were 2 sites that were missing - Stanford University and University of Miami. So I downloaded those sites from the [ABIDE-II webpage](http://fcon_1000.projects.nitrc.org/indi/abide/abide_II.html). I remaned them according to phenotype file: 

> ABIDEII-STANFORD -> ABIDEII-SU_2   
> ABIDEII-UM -> ABIDEII-U_MIA_1 

And converted them to BIDS format as follows:

## ABIDEII-UM 

### Anatomical File

```bash
python2 pathmatcher.py -ri "ABIDEII-UM/(\d{5})/session_(\d)/anat_(\d)/anat.nii.gz" -ro "ABIDEII-UM-BIDS/sub-\1/anat/sub-\1_T1w.nii.gz" -i ../ABIDE2RawDataBIDS/ -o ../ABIDE2RawDataBIDS/ -c 
```
  
### Functional Files
```bash
python2 pathmatcher.py -ri "ABIDEII-UM/(\d{5})/session_(\d)/rest_(\d)/rest.nii.gz" -ro "ABIDEII-UM-BIDS/sub-\1/func/sub-\1_task-rest_run-\3_bold.nii.gz" -i ../ABIDE2RawDataBIDS/ -o ../ABIDE2RawDataBIDS/ -c 
```
  

## ABIDEII-STANFORD 
### Anatommical Files
```bash
python2 pathmatcher.py -ri "ABIDEII-STANFORD/(\d{5})/session_(\d)/anat_(\d)/anat.nii.gz" -ro "ABIDEII-STANFORD-BIDS/sub-\1/anat/sub-\1_T1w.nii.gz" -i ../ABIDE2RawDataBIDS/ -o ../ABIDE2RawDataBIDS/ -c 
```
  
### Functional Files
```bash
python2 pathmatcher.py -ri "ABIDEII-STANFORD/(\d{5})/session_(\d)/rest_(\d)/rest.nii.gz" -ro "ABIDEII-STANFORD-BIDS/sub-\1/func/sub-\1_task-rest_run-\3_bold.nii.gz" -i ../ABIDE2RawDataBIDS/ -o ../ABIDE2RawDataBIDS/ -c 
```

Also, 2 other folders had different name in the phenotype file. So I renamed them as well.

> ABIDEII-ETHZ_1 -> ABIDEII-ETH_1   
> ABIDEII-ONRC_2 -> ABIDEII-OILH_2 


But I didn't use the above two sites in my analysis for the following reason.
The ABIDE-II data that I downloaded was in BIDS format but it didn have the metadata files (json files). I downloaded these files form [FCP-INDI](https://github.com/FCP-INDI/indi_bidsification/tree/master/ABIDE2/scan_jsons) and copied them into the respective sites.  
But it didn't have json files for the Stanford University and University of Miami sites so I decided not to include them at this time. Later I will create the metadata json files for these 2 sites and then use it.

## Issue 2


ABIDEII-SDSU_1/sub-28903/ses-1 contains readme_sub-28093_ses-1.txt. Note that the folder has subject ID 28903 whereas readme has got it wrong as 28093. So BIDS data grabber gets confused and while extracting the subject IDs:  
```bash
from bids.grabbids import BIDSLayout
layout = BIDSLayout(data_dir)    
subject_list = layout.get_subjects()
```
I extracted the subject ID 28093 which had no anatomical or functional file. So I had to rename the readme_sub-28093_ses-1.txt as readme_sub-28903_ses-1.txt

Also the following subject IDs had missing anatomical or functional file for Session=1 and Run=1:
```
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
```

So I discarded these subjects.



