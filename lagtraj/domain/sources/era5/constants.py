from scipy.constants import pi

# ECMWF CONSTANTS
rd = 287.06
rg = 9.80665
rv_over_rd_minus_one = 0.609133

# OTHER CONSTANTS NEED TO GO ELSEWHERE EVENTUALLY?
p_ref = 1.0e5
cp = 1004.0
rd_over_cp = rd / cp
p_ref_inv = 1.0 / p_ref
omega = 7.2921150e-5

longitude_tolerance = 0.001  # (about 100m)
