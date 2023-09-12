import numpy as np
from scipy import integrate
from T import T

def generateTestPositions(vel_expression, initial_position_3d, conditions):
    noise, delta_T, NOS, theta_degree, _, SRD, RDD,_,_ = conditions
    theta = theta_degree / 180 * np.pi

    r0_0 = initial_position_3d
    real_positions = np.zeros((NOS, 3))
    real_positions[0,:] = initial_position_3d
    M_p = (SRD + RDD) / (SRD + r0_0[2])

    x_proj = np.zeros(NOS)
    z_proj = np.zeros(NOS)

    x_proj[0] = M_p * r0_0[0] + np.random.randn() * noise / 2
    z_proj[0] = M_p * r0_0[2] + np.random.randn() * noise / 2

    for k in range(1, NOS):
        r0_k = r0_0 + np.array([
            integrate.quad(lambda t: vel_expression(t)[i], 0, delta_T * k)[0]
            for i in range(len(initial_position_3d))
        ])
        real_positions[k,:] = r0_k

        r_now = T(r0_k, theta * k)
        M_p = (SRD + RDD) / (SRD + r_now[1])

        x_proj[k] = M_p * r_now[0] + np.random.randn() * noise / 2
        z_proj[k] = M_p * r_now[2] + np.random.randn() * noise / 2

    xz_proj = np.column_stack((x_proj, z_proj))

    return xz_proj, real_positions

# Filename: particle_motion.py

import numpy as np
initial_position_A = [0.0, 0, 1]  # [x, y, z]
initial_position_B = [0.04, 0.01, 0]

# def velocity_A(t):
#     return lambda t: [np.sin(t), np.cos(t), 1]

# def velocity_B(t):
#     return lambda t: [np.sin(t + np.pi), np.cos(t + np.pi), -1]

def get_velocity_function(particle_id):
    if particle_id == 0:

        return lambda t: [np.sin(t), np.cos(t), 0.1]
    elif particle_id == 1:
        return lambda t: [0,0,0]
        return lambda t: [0.1, np.sin(t), np.cos(t)]
    
def get_initial_position(particle_id):
    switcher = {
        0: initial_position_A,
        1: initial_position_B
    }
    pos = switcher.get(particle_id)
    if pos:
        return pos
    else:
        raise ValueError("Invalid particle_id. Choose 1 for A and 2 for B.")

# Example usage:
if __name__ == "__main__":
    t = 5  # Some time value
    particle_id = 1  # Choose 1 for A and 2 for B

    vel_func = get_velocity_function(particle_id)
    velocity = vel_func(t)
    print(f"Velocity of Particle {particle_id} at time t:", velocity)
