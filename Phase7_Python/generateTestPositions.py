import numpy as np
from scipy import integrate
from T import T


# 

def generateTestPositions(vel_expression, initial_position_3d, conditions):
    noise, delta_T, NOS, theta_degree, _, SRD, RDD,_,_ = conditions
    theta = theta_degree / 180 * np.pi

    r0_0 = initial_position_3d
    real_positions = np.zeros((NOS, 3))
    real_positions[0,:] = initial_position_3d
    M_p = (SRD + RDD) / (SRD + r0_0[1])

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
# Modified initial positions for particles A and B
<<<<<<<< HEAD:Phase4/Variable_Speeds_Python/generateTestPositions.py
initial_position_A = [0.2, 0.2, 0]  # [x, y, z] (center of the semi-circle for particle A)
initial_position_B = [10, -15, 0]  # [x, y, z] (center of the semi-circle for particle B)
========
initial_position_A = [0, 0, 4]  # [x, y, z] (center of the semi-circle for particle A)
initial_position_B = [3, 5,-3]  # [x, y, z] (center of the semi-circle for particle B)
>>>>>>>> 337d6ad9fad7a8066a205f6769d529c5165f95ad:Phase7_Python/generateTestPositions.py
radius = 0.1  # Radius of the semi-circle

def get_velocity_function(particle_id):
    if particle_id == 0:
        # Particle A moves along a semi-circle going upwards in the xz-plane
        return lambda t: [radius * np.cos(np.pi * t), 0, radius * np.sin(np.pi * t)]
    elif particle_id == 1:
        # Particle B moves along a semi-circle going downwards in the xz-plane
        return lambda t: [-radius * np.cos(np.pi * t), 0, -radius * np.sin(np.pi * t)]
    else:
        raise ValueError("Invalid particle_id.")

def get_initial_position(particle_id):
    switcher = {
        0: initial_position_A,
        1: initial_position_B
    }
    pos = switcher.get(particle_id)
    if pos:
        return pos
    else:
        raise ValueError("Invalid particle_id. Choose 0 for A and 1 for B.")


# Example usage:
if __name__ == "__main__":
    t = 5  # Some time value
    particle_id = 1  # Choose 1 for A and 2 for B

    vel_func = get_velocity_function(particle_id)
    velocity = vel_func(t)
    print(f"Velocity of Particle {particle_id} at time t:", velocity)

