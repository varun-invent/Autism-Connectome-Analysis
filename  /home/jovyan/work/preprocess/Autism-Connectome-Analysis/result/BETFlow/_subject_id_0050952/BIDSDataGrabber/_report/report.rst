Node: BIDSDataGrabber (utility)
===============================

 Hierarchy : BETFlow.BIDSDataGrabber
 Exec ID : BIDSDataGrabber.a0

Original Inputs
---------------

* data_dir : /home/jovyan/work/preprocess/data/ABIDE-BIDS/NYU/
* function_str : def get_nifti_filenames(subject_id,data_dir):
#     Remember that all the necesary imports need to be INSIDE the function for the Function Interface to work!
    from bids.grabbids import BIDSLayout

    layout = BIDSLayout(data_dir)

    bold = [f.filename for f in layout.get(subject=subject_id, type='bold', extensions=['nii', 'nii.gz'])]

    return func_file_path

* ignore_exception : False
* subject_id : 0050952

