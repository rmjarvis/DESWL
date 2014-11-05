# This test measures the consistency of xi in SNR bins.
# It can be useful to assess the shear calibration errors.
# Basic workflow is listed below.
# To change parameters, including switching between ngmix and im3shape catalogs, see
# the config file test_xi_vs_snr.yaml.

# TODO: implement noise bias calibration/ responsivity corrections
# TODO: implement cuts

# get the weights
python test_xi_vs_snr.py  -c test_xi_vs_snr.yaml -a get_weights -n 394

# run the measurement
python test_xi_vs_snr.py  -c test_xi_vs_snr.yaml -a get_xi_vs_snr -n 394

# plot the result
python test_xi_vs_snr.py  -c test_xi_vs_snr.yaml -a plot_xi_vs_snr -n 394
