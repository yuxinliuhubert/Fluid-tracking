import numpy as np

def proj2r0_acc_combination(xz_proj, theta, SOD, ODD, dt):
    NOS = len(xz_proj)
    SDD = SOD + ODD
    row_number_A = round(2 * NOS + 2 * (np.math.factorial(NOS) / (np.math.factorial(NOS - 2) * 2)), 0)
    col_number_A = round(1 + 2 * NOS, 0)
    A = np.zeros((int(row_number_A), int(col_number_A))
    b = np.zeros((int(row_number_A), 1))

    # This loop is for constructing the equations arising from magnification alone
    for j in range(NOS):
        xi_j = xz_proj[j, 0]
        zi_j = xz_proj[j, 1]
        A[2 * j, 0] = 1
        A[2 * j, 2 * j + 1] = -zi_j / SDD
        b[2 * j] = zi_j * SOD / SDD
        A[2 * j + 1, 2 * j] = -1
        A[2 * j + 1, 2 * j + 1] = xi_j / SDD
        b[2 * j + 1] = -xi_j * SOD / SDD

    IoR = 2 * NOS + 1
    for k in range(2, NOS):
        trans_count = k - 1
        A[IoR:IoR + (trans_count * 2), 2 * k:2 * k + 2] = np.tile([[-1, 0], [0, -1]], (trans_count, 1))

        for i in range(trans_count):
            delta_theta = theta * i
            T1 = (k - 1) * dt - i * dt
            theta_prime = theta * (k - i - 1)
            A[IoR + 2 * i:IoR + 2 * i + 2, 2 * trans_count - 2 * i:2 * trans_count - 2 * i + 2] = np.array([
                [np.cos(delta_theta), np.sin(delta_theta)],
                [-np.sin(delta_theta), np.cos(delta_theta)]
            ])
            IoR += 2 * (k - 1)

    A = np.hstack([A, np.zeros(int(row_number_A), 6)])
    new_col_num = A.shape[1]
    u_ind = new_col_num - 5
    v_ind = u_ind + 1
    w_ind = u_ind + 2
    ax_ind = u_ind + 3
    ay_ind = u_ind + 4
    az_ind = u_ind + 5

    for j in range(NOS):
        A[2 * j, w_ind] = dt * (j)
        A[2 * j, az_ind] = 0.5 * (dt * (j) ** 2)

    IoR = 2 * NOS + 1

    for k in range(2, NOS):
        trans_count = k - 1
        T2 = dt * (k - 1)

        for i in range(trans_count):
            delta_theta = theta * i
            T1 = T2 - i * dt
            theta_prime = theta * (k - i - 1)
            A[IoR, u_ind] = (np.cos(theta_prime) * np.cos(delta_theta) - np.sin(theta_prime) * np.sin(delta_theta)) * (T2 - T1)
            A[IoR + 1, u_ind] = (-np.cos(theta_prime) * np.sin(delta_theta) - np.sin(theta_prime) * np.cos(delta_theta)) * (T2 - T1)
            A[IoR, v_ind] = (np.sin(theta_prime) * np.cos(delta_theta) + np.cos(theta_prime) * np.sin(delta_theta)) * (T2 - T1)
            A[IoR + 1, v_ind] = (-np.sin(theta_prime) * np.sin(delta_theta) + np.cos(theta_prime) * np.cos(delta_theta)) * (T2 - T1)
            A[IoR, ax_ind] = 0.5 * (np.cos(theta_prime) * np.cos(delta_theta) - np.sin(theta_prime) * np.sin(delta_theta)) * (T2 ** 2 - T1 ** 2)
            A[IoR + 1, ax_ind] = 0.5 * (-np.cos(theta_prime) * np.sin(delta_theta) - np.sin(theta_prime) * np.cos(delta_theta)) * (T2 ** 2 - T1 ** 2)
            A[IoR, ay_ind] = 0.5 * (np.sin(theta_prime) * np.cos(delta_theta) + np.cos(theta_prime) * np.sin(delta_theta)) * (T2 ** 2 - T1 ** 2)
            A[IoR + 1, ay_ind] = 0.5 * (-np.sin(theta_prime) * np.sin(delta_theta) + np.cos(theta_prime) * np.cos(delta_theta)) * (T2 ** 2 - T1 ** 2)
            IoR += 2

    x = np.linalg.lstsq(A, b, rcond=None)[0]
    r0 = [x[2], x[3], x[1], x[u_ind], x[v_ind], x[w_ind], x[ax_ind], x[ay_ind], x[az_ind]]
    return r0
