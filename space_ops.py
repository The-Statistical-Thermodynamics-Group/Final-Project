import numpy as np


def rotate_vector(vector, angle):
        z = complex(*vector) * np.exp(complex(0, angle))
        result = np.array([z.real, z.imag])
        return result


def get_norm(vector):
    return sum([x**2 for x in vector])**0.5
