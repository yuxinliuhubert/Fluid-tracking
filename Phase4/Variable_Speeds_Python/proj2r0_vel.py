import numpy as np


def proj2r0_acc_combination(xz_proj, theta, SOD, ODD, dt):
    NOS = xz_proj.shape[0]
    SDD = SOD + ODD
    row_number_A = int(2 * NOS + 2 * (np.math.factorial(NOS) / (np.math.factorial(NOS - 2) * 2)))
    col_number_A = int(1 + 2 * NOS)
    A = np.zeros((row_number_A, col_number_A))
    b = np.zeros((row_number_A, 1))

    for j in range(NOS):
        xi_j, zi_j = xz_proj[j]
        A[2 * j, 0] = 1
        A[2 * j, 2 * j +2] = -zi_j / SDD
        b[2 * j] = zi_j * SOD / SDD
        A[2 * j + 1, 2 * j+1] = -1
        A[2 * j + 1, 2 * j + 2] = xi_j / SDD
        b[2 * j + 1] = -xi_j * SOD / SDD

    IoR = 2 * NOS

    for k in range(1, NOS):
        trans_count = k
        A[IoR:IoR + 2 * trans_count, 2 * k+1:2 * k + 3] = np.tile(np.array([[-1, 0], [0, -1]]), (trans_count, 1))

        for i in range(trans_count):
            delta_theta = theta * (i + 1)
            A[IoR + 2 * i:IoR + 2 * i + 2, 2 * (k - 1) - 2 * i+1:2 * (k - 1) - 2 * i + 2+1] = np.array([[np.cos(delta_theta), np.sin(delta_theta)], [-np.sin(delta_theta), np.cos(delta_theta)]])
        IoR += 2 * trans_count

    A = np.hstack([A, np.zeros((A.shape[0], 3))])
    new_col_num = A.shape[1]
    u_ind, v_ind, w_ind = new_col_num - 3, new_col_num - 2, new_col_num - 1

    for j in range(NOS):
        A[2 * j, w_ind] = dt * j
        

    IoR = 2 * NOS

    for k in range(1, NOS):
        trans_count = k
        T2 = dt * k

        for i in range(trans_count):
            delta_theta = theta * (i + 1)
            T1 = T2 - (i+1) * dt
            theta_prime = theta * (k - i - 1)
            A[IoR, u_ind] = (np.cos(theta_prime) * np.cos(delta_theta) - np.sin(theta_prime) * np.sin(delta_theta)) * (T2 - T1)
            A[IoR + 1, u_ind] = (-np.cos(theta_prime) * np.sin(delta_theta) - np.sin(theta_prime) * np.cos(delta_theta)) * (T2 - T1)
            A[IoR, v_ind] = (np.sin(theta_prime) * np.cos(delta_theta) + np.cos(theta_prime) * np.sin(delta_theta)) * (T2 - T1)
            A[IoR + 1, v_ind] = (-np.sin(theta_prime) * np.sin(delta_theta) + np.cos(theta_prime) * np.cos(delta_theta)) * (T2 - T1)
            IoR += 2

    x, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)
    r0 = [x[1], x[2], x[0], x[u_ind], x[v_ind], x[w_ind], 0,0,0]
    return r0
