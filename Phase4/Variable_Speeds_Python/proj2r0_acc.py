import numpy as np

def proj2r0_acc(proj, theta, SRD, RDD, delta_T):
    NOS = len(proj)
    SDD = SRD + RDD
    row_number_A = 2 * NOS + 2 * (NOS - 1)
    col_number_A = 1 + 2 * NOS
    A = np.zeros((row_number_A, col_number_A))
    b = np.zeros(row_number_A)
    
    for j in range(NOS):
        xi_j, zi_j = proj[j]
        A[2*j, 0] = 1
        A[2*j, 2*j+2] = -zi_j / SDD
        b[2*j] = zi_j * SRD / SDD
        
        A[2*j+1, 2*j+1] = -1
        A[2*j+1, 2*j+2] = xi_j / SDD
        b[2*j+1] = -xi_j * SRD / SDD
    
    x = 2 * NOS
    for k in range(1, NOS):
        A[x:x+2, 2*k+1:2*k+3] = [[-1, 0], [0, -1]]
        A[x:x+2, 1:3] = [[np.cos(theta*(k)), np.sin(theta*(k))], [-np.sin(theta*(k)), np.cos(theta*(k))]]
        x += 2
        
    A = np.pad(A, ((0, 0), (0, 6)), 'constant')
    new_col_num = A.shape[1]
    
    for j in range(NOS):
        A[2*j, new_col_num-4] = delta_T * (j)
        A[2*j, new_col_num-1] = 0.5 * (delta_T * (j))**2
        
    IoR = 2 * NOS
    for k in range(1, NOS):
        A[IoR, new_col_num-6] = np.cos(theta * k) * delta_T * k
        A[IoR+1, new_col_num-6] = -np.sin(theta * k) * delta_T * k
        
        A[IoR, new_col_num-5] = np.sin(theta * k) * delta_T * k
        A[IoR+1, new_col_num-5] = np.cos(theta * k) * delta_T * k
        
        A[IoR, new_col_num-3] = 0.5 * np.cos(theta * k) * (delta_T * k)**2
        A[IoR+1, new_col_num-3] = -np.sin(theta * k) * 0.5 * (delta_T * k)**2
        
        A[IoR, new_col_num-2] = 0.5 * np.sin(theta * k) * (delta_T * k)**2
        A[IoR+1, new_col_num-2] = 0.5 * np.cos(theta * k) * (delta_T * k)**2
        IoR += 2
        
    x, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)
    r0 = [x[1], x[2], x[0], x[new_col_num-6], x[new_col_num-5], x[new_col_num-4], x[new_col_num-3], x[new_col_num-2], x[new_col_num-1]]
    return r0

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

    A = np.hstack([A, np.zeros((A.shape[0], 6))])
    new_col_num = A.shape[1]
    u_ind, v_ind, w_ind, ax_ind, ay_ind, az_ind = new_col_num - 6, new_col_num - 5, new_col_num - 4, new_col_num - 3, new_col_num - 2, new_col_num - 1

    for j in range(NOS):
        A[2 * j, w_ind] = dt * j
        A[2 * j, az_ind] = 0.5 * (dt * j)**2

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
            A[IoR, ax_ind] = 0.5 * (np.cos(theta_prime) * np.cos(delta_theta) - np.sin(theta_prime) * np.sin(delta_theta)) * (T2**2 - T1**2)
            A[IoR + 1, ax_ind] = 0.5 * (-np.cos(theta_prime) * np.sin(delta_theta) - np.sin(theta_prime) * np.cos(delta_theta)) * (T2**2 - T1**2)
            A[IoR, ay_ind] = 0.5 * (np.sin(theta_prime) * np.cos(delta_theta) + np.cos(theta_prime) * np.sin(delta_theta)) * (T2**2 - T1**2)
            A[IoR + 1, ay_ind] = 0.5 * (-np.sin(theta_prime) * np.sin(delta_theta) + np.cos(theta_prime) * np.cos(delta_theta)) * (T2**2 - T1**2)
            IoR += 2

    x, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)
    r0 = [x[1], x[2], x[0], x[u_ind], x[v_ind], x[w_ind], x[ax_ind], x[ay_ind], x[az_ind]]
    return r0

def proj2r0_vel(xz_proj, theta, SOD, ODD, dt):
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

def proj2r0_stationary(xz_proj, theta, SOD, ODD, dt):
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

 
    x, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)
    r0 = [x[1], x[2], x[0], 0,0,0, 0,0,0]
    return r0