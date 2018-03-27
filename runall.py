import postprocessingfc as ppfc
import itertools

itr = (list(itertools.product([0, 1], repeat=3)))

itr = [(1,0,1,1)]
number_of_subjects = -1

for motion_param_regression, global_signal_regression, band_pass_filtering, smoothing in itr:
    combination = 'motionRegress' + str(int(motion_param_regression)) + 'filt' + \
              str(int(band_pass_filtering)) + 'global' + str(int(global_signal_regression)) + \
              'smoothing' + str(int(smoothing))

    print("Combination: ",combination)
    functional_connectivity_directory =  combination
    print(motion_param_regression, global_signal_regression, band_pass_filtering, smoothing)
    ppfc.main(motion_param_regression, global_signal_regression, band_pass_filtering, smoothing, number_of_subjects,functional_connectivity_directory)
