import numpy as np
from scipy.constants import pi

from ... import njit


@njit
def steffen_3d(
    v_in,
    z_in,
    z_out,
    z_min_surface,
    z_max_surface,
    lower_extrapolation_with_gradient=False,
    upper_extrapolation_with_gradient=False,
):
    """
    Performs Steffen interpolation on each individual column.
    Steffen, M. (1990). A simple method for monotonic interpolation
    in one dimension. Astronomy and Astrophysics, 239, 443.

    `z_min_surface` and `z_max_surface` define the minimum and maximum height
    outside of which the interpolated value will be set to nan (useful for not
    interpolating pressure into solid ground for example).
    `lower_extrapolation_with_gradient` and `upper_extrapolation_with_gradient`
    enables the use of one-sided gradient when extrapolating outside the range
    of `z_in`, otherwise the limit value of `v_in` is used by default
    """
    assert v_in.shape[1:] == z_min_surface.shape
    k_max = v_in.shape[0]
    j_max = v_in.shape[1]
    i_max = v_in.shape[2]
    k_max_output = z_out.shape[0]
    k_max_minus = k_max - 1
    linear_slope = np.empty((k_max))
    v_out = np.empty((k_max_output, j_max, i_max))
    for i in range(i_max):
        for j in range(j_max):
            # first point
            delta_lower = z_in[1, j, i] - z_in[0, j, i]
            delta_upper = z_in[2, j, i] - z_in[1, j, i]
            if delta_lower < 0:
                raise Exception("Non-montonic increase in z_in")
            if delta_upper < 0:
                raise Exception("Non-montonic increase in z_in")
            slope_lower = (v_in[1, j, i] - v_in[0, j, i]) / delta_lower
            slope_upper = (v_in[2, j, i] - v_in[1, j, i]) / delta_upper

            weighted_slope = slope_lower * (
                1 + delta_lower / (delta_lower + delta_upper)
            ) - slope_upper * delta_lower / (delta_lower + delta_upper)
            if weighted_slope * slope_lower <= 0.0:
                linear_slope[0] = 0.0
            elif np.abs(weighted_slope) > 2 * np.abs(slope_lower):
                linear_slope[0] = 2.0 * slope_lower
            else:
                linear_slope[0] = weighted_slope

            # intermediate points
            for k in range(1, k_max_minus):
                delta_lower = z_in[k, j, i] - z_in[k - 1, j, i]
                delta_upper = z_in[k + 1, j, i] - z_in[k, j, i]
                slope_lower = (v_in[k, j, i] - v_in[k - 1, j, i]) / delta_lower
                slope_upper = (v_in[k + 1, j, i] - v_in[k, j, i]) / delta_upper
                weighted_slope = (
                    slope_lower * delta_upper + slope_upper * delta_lower
                ) / (delta_lower + delta_upper)

                if slope_lower * slope_upper <= 0.0:
                    linear_slope[k] = 0.0
                elif np.abs(weighted_slope) > 2.0 * np.abs(slope_lower):
                    linear_slope[k] = np.copysign(2.0, slope_lower) * min(
                        np.abs(slope_lower), np.abs(slope_upper)
                    )
                elif np.abs(weighted_slope) > 2.0 * np.abs(slope_upper):
                    linear_slope[k] = np.copysign(2.0, slope_lower) * min(
                        np.abs(slope_lower), np.abs(slope_upper)
                    )
                else:
                    linear_slope[k] = weighted_slope

            # last point
            delta_lower = z_in[k_max_minus - 1, j, i] - z_in[k_max_minus - 2, j, i]
            delta_upper = z_in[k_max_minus, j, i] - z_in[k_max_minus - 1, j, i]
            slope_lower = (
                v_in[k_max_minus - 1, j, i] - v_in[k_max_minus - 2, j, i]
            ) / delta_lower
            slope_upper = (
                v_in[k_max_minus, j, i] - v_in[k_max_minus - 1, j, i]
            ) / delta_upper
            weighted_slope = slope_upper * (
                1 + delta_upper / (delta_upper + delta_lower)
            ) - slope_lower * delta_upper / (delta_upper + delta_lower)
            if weighted_slope * slope_upper <= 0.0:
                linear_slope[k_max_minus] = 0.0
            elif np.abs(weighted_slope) > 2.0 * np.abs(slope_upper):
                linear_slope[k_max_minus] = 2.0 * slope_upper
            else:
                linear_slope[k_max_minus] = weighted_slope

            # loop over output points
            k_temp = 0
            for k_out in range(k_max_output):
                while (k_temp < k_max) and (z_in[k_temp, j, i] < z_out[k_out]):
                    k_temp = k_temp + 1
                if 0 < k_temp < k_max:
                    k_high = k_temp
                    k_low = k_high - 1
                    delta = z_in[k_high, j, i] - z_in[k_low, j, i]
                    slope = (v_in[k_high, j, i] - v_in[k_low, j, i]) / delta
                    a = (linear_slope[k_low] + linear_slope[k_high] - 2 * slope) / (
                        delta * delta
                    )
                    b = (
                        3 * slope - 2 * linear_slope[k_low] - linear_slope[k_high]
                    ) / delta
                    c = linear_slope[k_low]
                    d = v_in[k_low, j, i]
                    t_1 = z_out[k_out] - z_in[k_low, j, i]
                    t_2 = t_1 * t_1
                    t_3 = t_2 * t_1
                    v_out[k_out, j, i] = a * t_3 + b * t_2 + c * t_1 + d
                elif (k_temp == 0) and (z_out[k_out] >= z_min_surface[j, i]):
                    if lower_extrapolation_with_gradient:
                        v_out[k_out, j, i] = v_in[0, j, i] + linear_slope[0] * (
                            z_out[k_out] - z_in[0, j, i]
                        )
                    else:
                        v_out[k_out, j, i] = v_in[0, j, i]
                elif (k_temp == k_max) and (z_out[k_out] <= z_max_surface[j, i]):
                    if upper_extrapolation_with_gradient:
                        v_out[k_out, j, i] = v_in[k_max - 1, j, i] + linear_slope[
                            k_max - 1
                        ] * (z_out[k_out] - z_in[k_max - 1, j, i])
                    else:
                        v_out[k_out, j, i] = v_in[k_max - 1, j, i]
                else:
                    v_out[k_out, j, i] = np.nan
    return v_out


# This can probably be replaced by a generic Steffen interpolation function
# No extrapolation is performed
# Time axis present, but not lat, lon
@njit
def steffen_1d_no_ep_time(
    input_data, input_levels, output_level_array,
):
    """ Performs Steffen interpolation on one individual column for
    each time step.
    Steffen, M. (1990). A simple method for monotonic interpolation in
    one dimension. Astronomy and Astrophysics, 239, 443. """
    t_max = input_data.shape[0]
    k_max = input_data.shape[1]
    k_max_output = output_level_array.shape[0]
    k_max_minus = k_max - 1
    linear_slope = np.empty((k_max))
    output_data = np.empty((t_max, k_max_output))
    # first point
    delta_lower = input_levels[1] - input_levels[0]
    delta_upper = input_levels[2] - input_levels[1]
    if delta_lower < 0:
        raise Exception("Non-montonic increase in input_levels")
    if delta_upper < 0:
        raise Exception("Non-montonic increase in input_levels")
    for time_index in range(t_max):
        slope_lower = (
            input_data[time_index, 1] - input_data[time_index, 0]
        ) / delta_lower
        slope_upper = (
            input_data[time_index, 2] - input_data[time_index, 1]
        ) / delta_upper
        weighted_slope = slope_lower * (
            1 + delta_lower / (delta_lower + delta_upper)
        ) - slope_upper * delta_lower / (delta_lower + delta_upper)
        if weighted_slope * slope_lower <= 0.0:
            linear_slope[0] = 0.0
        elif np.abs(weighted_slope) > 2 * np.abs(slope_lower):
            linear_slope[0] = 2.0 * slope_lower
        else:
            linear_slope[0] = weighted_slope

        # intermediate points
        for k in range(1, k_max_minus):
            delta_lower = input_levels[k] - input_levels[k - 1]
            delta_upper = input_levels[k + 1] - input_levels[k]
            slope_lower = (
                input_data[time_index, k] - input_data[time_index, k - 1]
            ) / delta_lower
            slope_upper = (
                input_data[time_index, k + 1] - input_data[time_index, k]
            ) / delta_upper
            weighted_slope = (slope_lower * delta_upper + slope_upper * delta_lower) / (
                delta_lower + delta_upper
            )

            if slope_lower * slope_upper <= 0.0:
                linear_slope[k] = 0.0
            elif np.abs(weighted_slope) > 2.0 * np.abs(slope_lower):
                linear_slope[k] = np.copysign(2.0, slope_lower) * min(
                    np.abs(slope_lower), np.abs(slope_upper)
                )
            elif np.abs(weighted_slope) > 2.0 * np.abs(slope_upper):
                linear_slope[k] = np.copysign(2.0, slope_lower) * min(
                    np.abs(slope_lower), np.abs(slope_upper)
                )
            else:
                linear_slope[k] = weighted_slope

        # last point
        delta_lower = input_levels[k_max_minus - 1] - input_levels[k_max_minus - 2]
        delta_upper = input_levels[k_max_minus] - input_levels[k_max_minus - 1]
        slope_lower = (
            input_data[time_index, k_max_minus - 1]
            - input_data[time_index, k_max_minus - 2]
        ) / delta_lower
        slope_upper = (
            input_data[time_index, k_max_minus]
            - input_data[time_index, k_max_minus - 1]
        ) / delta_upper
        weighted_slope = slope_upper * (
            1 + delta_upper / (delta_upper + delta_lower)
        ) - slope_lower * delta_upper / (delta_upper + delta_lower)
        if weighted_slope * slope_upper <= 0.0:
            linear_slope[k_max_minus] = 0.0
        elif np.abs(weighted_slope) > 2.0 * np.abs(slope_upper):
            linear_slope[k_max_minus] = 2.0 * slope_upper
        else:
            linear_slope[k_max_minus] = weighted_slope

        # loop over output points
        k_temp = 0
        for k_out in range(k_max_output):
            while (k_temp < k_max) and (
                input_levels[k_temp] < output_level_array[k_out]
            ):
                k_temp = k_temp + 1
            if 0 < k_temp < k_max:
                k_high = k_temp
                k_low = k_high - 1
                delta = input_levels[k_high] - input_levels[k_low]
                slope = (
                    input_data[time_index, k_high] - input_data[time_index, k_low]
                ) / delta
                a = (linear_slope[k_low] + linear_slope[k_high] - 2 * slope) / (
                    delta * delta
                )
                b = (3 * slope - 2 * linear_slope[k_low] - linear_slope[k_high]) / delta
                c = linear_slope[k_low]
                d = input_data[time_index, k_low]
                t_1 = output_level_array[k_out] - input_levels[k_low]
                t_2 = t_1 * t_1
                t_3 = t_2 * t_1
                output_data[time_index, k_out] = a * t_3 + b * t_2 + c * t_1 + d
            # Allow for small deviations in upper/lower levels
            elif (k_temp == 0) and (
                abs(output_level_array[k_out] - input_levels[k_temp]) < 1e-6
            ):
                output_data[time_index, k_out] = input_data[time_index, 0]
            elif (k_temp == k_max) and (
                abs(output_level_array[k_out] - input_levels[k_temp]) < 1e-6
            ):
                output_data[time_index, k_out] = input_data[time_index, k_max]
            else:
                output_data[time_index, k_out] = np.nan
    return output_data


def central_estimate(a_in):
    """take a one-sided difference at the edges, and a central difference elsewhere"""
    return np.concatenate(([a_in[0]], 0.5 * (a_in[1:-1] + a_in[2:]), [a_in[-1]]))


def cos_transition(absolute_input, transition_start, transition_end):
    """function that smoothly transitions from 1 to 0 using a
    cosine-shaped transition between start and end"""
    normalised_input = (absolute_input - transition_start) / (
        transition_end - transition_start
    )
    weight_factor = 1.0 * (normalised_input < 0.0) + (
        0.5 + 0.5 * np.cos(normalised_input * pi)
    ) * (1.0 - (normalised_input < 0.0) - (normalised_input > 1.0))
    return weight_factor
