import numpy as np
import T

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
