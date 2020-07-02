# Optional numba dependency
try:
    from numba import njit

    print("Running with numba")
except ImportError:
    def njit(numba_function):
        """Dummy numba function"""
        return numba_function

    print("Running without numba")


@njit
def steffen_3d(
    input_data,
    input_levels,
    output_level_array,
    lower_extrapolation_surface,
    lower_extrapolation_with_gradient=False,
):
    """
    Performs Steffen interpolation on each individual column.
    Steffen, M. (1990). A simple method for monotonic interpolation
    in one dimension. Astronomy and Astrophysics, 239, 443.
    """
    k_max = input_data.shape[0]
    j_max = input_data.shape[1]
    i_max = input_data.shape[2]
    k_max_output = output_level_array.shape[0]
    k_max_minus = k_max - 1
    linear_slope = np.empty((k_max, j_max, i_max))
    output_data = np.empty((k_max_output, j_max, i_max))
    for i in range(i_max):
        for j in range(j_max):
            # first point
            delta_lower = input_levels[1, j, i] - input_levels[0, j, i]
            delta_upper = input_levels[2, j, i] - input_levels[1, j, i]
            if delta_lower < 0:
                raise Exception("Non-montonic increase in input_levels")
            if delta_upper < 0:
                raise Exception("Non-montonic increase in input_levels")
            slope_lower = (input_data[1, j, i] - input_data[0, j, i]) / delta_lower
            slope_upper = (input_data[2, j, i] - input_data[1, j, i]) / delta_upper

            weighted_slope = slope_lower * (
                1 + delta_lower / (delta_lower + delta_upper)
            ) - slope_upper * delta_lower / (delta_lower + delta_upper)
            if weighted_slope * slope_lower <= 0.0:
                linear_slope[0, j, i] = 0.0
            elif np.abs(weighted_slope) > 2 * np.abs(slope_lower):
                linear_slope[0, j, i] = 2.0 * slope_lower
            else:
                linear_slope[0, j, i] = weighted_slope

            # intermediate points
            for k in range(1, k_max_minus):
                delta_lower = input_levels[k, j, i] - input_levels[k - 1, j, i]
                delta_upper = input_levels[k + 1, j, i] - input_levels[k, j, i]
                slope_lower = (
                    input_data[k, j, i] - input_data[k - 1, j, i]
                ) / delta_lower
                slope_upper = (
                    input_data[k + 1, j, i] - input_data[k, j, i]
                ) / delta_upper
                weighted_slope = (
                    slope_lower * delta_upper + slope_upper * delta_lower
                ) / (delta_lower + delta_upper)

                if slope_lower * slope_upper <= 0.0:
                    linear_slope[k, j, i] = 0.0
                elif np.abs(weighted_slope) > 2.0 * np.abs(slope_lower):
                    linear_slope[k, j, i] = np.copysign(2.0, slope_lower) * min(
                        np.abs(slope_lower), np.abs(slope_upper)
                    )
                elif np.abs(weighted_slope) > 2.0 * np.abs(slope_upper):
                    linear_slope[k, j, i] = np.copysign(2.0, slope_lower) * min(
                        np.abs(slope_lower), np.abs(slope_upper)
                    )
                else:
                    linear_slope[k, j, i] = weighted_slope

            # last point
            delta_lower = (
                input_levels[k_max_minus - 1, j, i]
                - input_levels[k_max_minus - 2, j, i]
            )
            delta_upper = (
                input_levels[k_max_minus, j, i] - input_levels[k_max_minus - 1, j, i]
            )
            slope_lower = (
                input_data[k_max_minus - 1, j, i] - input_data[k_max_minus - 2, j, i]
            ) / delta_lower
            slope_upper = (
                input_data[k_max_minus, j, i] - input_data[k_max_minus - 1, j, i]
            ) / delta_upper
            weighted_slope = slope_upper * (
                1 + delta_upper / (delta_upper + delta_lower)
            ) - slope_lower * delta_upper / (delta_upper + delta_lower)
            if weighted_slope * slope_upper <= 0.0:
                linear_slope[k_max_minus, j, i] = 0.0
            elif np.abs(weighted_slope) > 2.0 * np.abs(slope_upper):
                linear_slope[k_max_minus, j, i] = 2.0 * slope_upper
            else:
                linear_slope[k_max_minus, j, i] = weighted_slope

            # loop over output points
            k_temp = 0
            for k_out in range(k_max_output):
                while (k_temp < k_max) and (
                    input_levels[k_temp, j, i] < output_level_array[k_out]
                ):
                    k_temp = k_temp + 1
                if 0 < k_temp < k_max:
                    k_high = k_temp
                    k_low = k_high - 1
                    delta = input_levels[k_high, j, i] - input_levels[k_low, j, i]
                    slope = (input_data[k_high, j, i] - input_data[k_low, j, i]) / delta
                    a = (
                        linear_slope[k_low, j, i]
                        + linear_slope[k_high, j, i]
                        - 2 * slope
                    ) / (delta * delta)
                    b = (
                        3 * slope
                        - 2 * linear_slope[k_low, j, i]
                        - linear_slope[k_high, j, i]
                    ) / delta
                    c = linear_slope[k_low, j, i]
                    d = input_data[k_low, j, i]
                    t_1 = output_level_array[k_out] - input_levels[k_low, j, i]
                    t_2 = t_1 * t_1
                    t_3 = t_2 * t_1
                    output_data[k_out, j, i] = a * t_3 + b * t_2 + c * t_1 + d
                elif (k_temp == 0) and (
                    output_level_array[k_out] >= lower_extrapolation_surface[j, i]
                ):
                    if lower_extrapolation_with_gradient:
                        output_data[k_out, j, i] = input_data[0, j, i] + linear_slope[
                            0, j, i
                        ] * (output_level_array[k_out] - input_levels[0, j, i])
                    else:
                        output_data[k_out, j, i] = input_data[0, j, i]
                else:
                    output_data[k_out, j, i] = np.nan
    return output_data


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
