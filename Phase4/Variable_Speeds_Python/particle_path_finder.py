import numpy as np
from T import T
def track_particles(data_set, alpha):
    M, N = data_set.shape[0], data_set.shape[1] // 2  # Assuming each particle has 2 coordinates (x, z)
    tracked_particles = np.zeros((M, 2 * N))
    tracked_particles[0] = data_set[0]

    for m in range(1, M):
        for n in range(0, 2 * N, 2):
            particle = data_set[m, n:n+2]
            # rotated_particle = T(particle, alpha)
            previous_shot = tracked_particles[m-1].reshape(N, 2)
            closest_particle_index = find_closest_particle(particle, previous_shot)
            tracked_particles[m, 2*closest_particle_index:2*closest_particle_index+2] = particle

    return tracked_particles

def find_closest_particle(particle, shot):
    distances = np.linalg.norm(shot - particle, axis=1)
    return np.argmin(distances)

