import numpy as np
from scipy.constants import pi


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
