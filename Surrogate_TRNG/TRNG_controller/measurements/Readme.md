-> Bit_read_1M
    Contains the TRN bits collected from the surrogate under open-loop and closed-loop conditions
    The codes are available (ESP code to set the controller, and Python notebook script to automate the readings). The conditions under which the bits were collected are marked in the filenames. Data files are in NumPy format (.npz).

-> Hysteresis
    Hysteresis measured at different frequencies. Analysis codes are included.

-> Ktime
    Kramers time measured for 3 different bias conditions.

-> P_control_noise_variation
    Variation of noise power while keeping the probability set value fixed.

-> P_control_Pset_variation
    Variation of the probability set value while keeping the noise power fixed.

-> Probability
    Measurement of the corresponding P curves for the measured Kramers time (by varying the set
    pulse duration).

-> Probability_vs_bias
    Measurement of P curves by varying the bias values.

All the corresponding measured data are ploted in the controller_data.ipynb notebook in the project folder.