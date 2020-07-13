import numpy as np

try:
    from numba import njit
except ImportError:

    def njit(numba_function):
        """Dummy numba function"""
        return numba_function


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
