import postprocessingfc as ppfc
import itertools

itr = (list(itertools.product([0, 1], repeat=3)))
number_of_subjects = -1
for motion_param_regression, global_signal_regression, band_pass_filtering in itr:
    combination = 'fc_motionRegress' + str(int(motion_param_regression)) + 'filt' + str(int(band_pass_filtering)) + 'global' + str(int(global_signal_regression))
    print("Combination: ",combination)
    functional_connectivity_directory =  combination
    print(motion_param_regression, global_signal_regression, band_pass_filtering)
    ppfc.main(motion_param_regression, global_signal_regression, band_pass_filtering,number_of_subjects,functional_connectivity_directory)
