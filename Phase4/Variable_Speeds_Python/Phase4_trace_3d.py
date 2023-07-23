import numpy as np
from scipy import integrate
from scipy.stats import norm
from proj2r0_acc import proj2r0_acc
from T import T

def Phase4_trace_3d(initial_position_3d, conditions, vel_expression, xz_proj):
    _,delta_T, NOS, theta_degree, N, SRD, RDD,method,dataPiling = conditions

    theta = np.deg2rad(theta_degree)

    proj_used_index = 0
    
    NOS = int(conditions[2])  # Convert to int before using
    positions_predicted = np.zeros((NOS, 3))
    
    NOS_per_section = N
    prev_NOS_section = NOS_per_section
    print("NOS: ",NOS)
    print("NOS_per_Section: ",NOS_per_section) 

    if dataPiling == 'serial':
        while proj_used_index < NOS:
            alpha = -theta*(proj_used_index)
            temp = generateEstimatedPositions(alpha, proj_used_index, NOS_per_section,xz_proj,conditions)
            positions_predicted[proj_used_index : proj_used_index+NOS_per_section, :] = temp
            # print("normal called")

            proj_used_index += NOS_per_section 
            # proj_used_index += 1

            if abs(NOS - proj_used_index) < N:
                # NOS_per_section = NOS - proj_used_index + 1
                # alpha = -theta*(proj_used_index)
                # temp = generateEstimatedPositions(alpha, proj_used_index, NOS_per_section-1,xz_proj,conditions)
                # positions_predicted[proj_used_index : proj_used_index+NOS_per_section, :] = temp
                # proj_used_index += NOS_per_section

                prev_proj_index = proj_used_index
                proj_used_index = NOS - NOS_per_section
                alpha = -theta*(proj_used_index)
                last_positions = positions_predicted[proj_used_index : prev_proj_index, :]
                new_positions = generateEstimatedPositions(alpha, proj_used_index, NOS_per_section,xz_proj,conditions)
                print('prev',prev_proj_index)
                print('proj_index',proj_used_index)
                print('last positions',last_positions)
                print('new positions',new_positions)
                combined_positions = np.concatenate([(new_positions[:len(last_positions), :] + last_positions)/2, new_positions[len(last_positions):, :]], axis=0)

                positions_predicted[proj_used_index : proj_used_index+NOS_per_section, :] = combined_positions
                # print("retro called")
                proj_used_index += NOS_per_section + 1
                

    elif dataPiling == 'overlap':
        for i in range(int(NOS - N)):
            alpha = -theta * proj_used_index

            if proj_used_index == 0:
                positions_predicted[proj_used_index : N, :] = generateEstimatedPositions(alpha, proj_used_index, NOS_per_section,xz_proj,conditions)
            else:
                last_positions = positions_predicted[proj_used_index : proj_used_index + prev_NOS_section - 1, :]
                new_positions = generateEstimatedPositions(alpha, proj_used_index, NOS_per_section,xz_proj,conditions)
                new_positions = np.concatenate([(new_positions[:len(last_positions), :] + last_positions)/2, new_positions[len(last_positions):, :]], axis=0)
                positions_predicted[proj_used_index : proj_used_index+len(new_positions)-1, :] = new_positions
            
            proj_used_index += 1

    return positions_predicted

def generateEstimatedPositions(alpha, proj_used_index, N, xz_proj, conditions):
    _,delta_T, _, theta_degree, _, SRD, RDD,method,dataPiling = conditions
    theta = np.deg2rad(theta_degree)
    positions_predicted = np.zeros((N, 3))

    values_this_round = proj2r0_acc(xz_proj[proj_used_index : proj_used_index+N-2, :], theta, SRD, RDD, delta_T)

    if method == 'acceleration':
        x0, y0, z0, u, v, w, a_x, a_y, a_z = values_this_round
        position_rotated =  np.transpose(T([x0, y0, z0], alpha))
        x0, y0, z0 = position_rotated
        positions_predicted[0, :] = position_rotated
        velocity_rotated = np.transpose(T([u, v, w], alpha))
        u, v, w = velocity_rotated
        acc_rotated =  np.transpose(T([a_x, a_y, a_z], alpha))
        a_x, a_y, a_z = acc_rotated

        for j in range(1, N):
            time = delta_T * (j)
            positions_predicted[j, :] = [x0+u*time+0.5*a_x*time**2, y0+v*time+0.5*a_y*time**2, z0+w*time+0.5*a_z*time**2]

    elif method == 'linear':
        x0, y0, z0, u, v, w = values_this_round
        position_rotated = np.transpose(T([x0, y0, z0], alpha))
        x0, y0, z0 = position_rotated
        positions_predicted[0, :] = position_rotated
        velocity_rotated = np.transpose(T([u, v, w], alpha))
        u, v, w = velocity_rotated

        for j in range(1, N):
            time = delta_T * (j)
            positions_predicted[j, :] = [x0+u*time, y0+v*time, z0+w*time]

    return positions_predicted
