import numpy as np

def T(r1, alpha):
    rotation_matrix = np.array([
        [np.cos(-alpha), -np.sin(-alpha), 0],
        [np.sin(-alpha),  np.cos(-alpha), 0],
        [0,               0,              1]
    ])
    r2 = np.matmul(rotation_matrix, r1)
    return r2
