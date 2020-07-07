import numpy as np


# Optional numba dependency
try:
    from numba import njit

    print("Running with numba")
except ImportError:

    def njit(numba_function):
        """Dummy numba function"""
        return numba_function

    print("Running without numba")


def cos_transition(absolute_input, transition_start, transition_end):
    """function that smoothly transitions from 1 to 0 using a cosine-shaped
    transition between start and end"""
    normalised_input = (absolute_input - transition_start) / (
        transition_end - transition_start
    )
    weight_factor = 1.0 * (normalised_input < 0.0) + (
        0.5 + 0.5 * np.cos(normalised_input * pi)
    ) * (1.0 - (normalised_input < 0.0) - (normalised_input > 1.0))
    return weight_factor


@njit
def boundary_gradients(x_array, y_array, val_array):
    """Numba function to calculate gradients, dismissing filtered points"""
    len_temp = np.shape(val_array)[0]
    len_levels = np.shape(val_array)[1]
    len_lats = np.shape(val_array)[2]
    len_lons = np.shape(val_array)[3]
    x_gradient_array = np.empty((len_temp, len_levels))
    y_gradient_array = np.empty((len_temp, len_levels))
    dval_dx = np.empty((len_lats))
    dval_dy = np.empty((len_lons))
    for this_time in range(len_temp):
        for this_level in range(len_levels):
            for this_lat in range(len_lats):
                vals_at_lat = val_array[this_time, this_level, this_lat, :].flatten()
                x_at_lat = x_array[this_lat, :].flatten()
                vals_filtered = vals_at_lat[~np.isnan(vals_at_lat)]
                x_filtered = x_at_lat[~np.isnan(vals_at_lat)]
                dvals = vals_filtered[-1] - vals_filtered[0]
                dval_dx[this_lat] = dvals / (x_filtered[-1] - x_filtered[0])
            for this_lon in range(len_lons):
                vals_at_lon = val_array[this_time, this_level, :, this_lon].flatten()
                y_at_lat = y_array[:, this_lon].flatten()
                vals_filtered = vals_at_lon[~np.isnan(vals_at_lon)]
                y_filtered = y_at_lat[~np.isnan(vals_at_lon)]
                dvals = vals_filtered[-1] - vals_filtered[0]
                dval_dy[this_lon] = dvals / (y_filtered[-1] - y_filtered[0])
            x_gradient_array[this_time, this_level] = np.mean(
                dval_dx[~np.isnan(dval_dx)]
            )
            y_gradient_array[this_time, this_level] = np.mean(
                dval_dy[~np.isnan(dval_dy)]
            )
    return x_gradient_array, y_gradient_array
