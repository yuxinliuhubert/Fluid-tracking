import numpy as np
from T import T  # Assuming this function rotates a particle by an angle 'alpha'
import heapq
from Phase4_trace_3d import Phase4_trace_3d
# make sure we make enough defensive copies of the data
class ParticlePathFinder:


    def __init__(self, alpha, conditions, reconstruction_conditions) -> None:
        """
        Initializes a ParticlePathFinder object with the given parameters.

        Args:
        - alpha (float): the alpha value used in the particle path finding algorithm
        - conditions (list): a list of conditions used in the particle path finding algorithm
            - learning_rate_2D (float): the learning rate used in the 2D particle path finding algorithm
            - motion_randomness (float): the motion randomness used in the particle path finding algorithm
            - learning_rate_3D (float): the learning rate used in the 3D particle path finding algorithm
        - reconstruction_conditions (list): a list of conditions used in the particle path finding algorithm
            - NOS (int): the number of sections used in the particle path finding algorithm
            - NOS_per_section (int): the number of particles per section used in the particle path finding algorithm
        """

        self.alpha = alpha
        self.current_snapShotIndex = 0
        self.particle_id = 0
        self.particleData_2D = {}
        self.shotData = {}
        self.shotData_original = {}
        self.learning_rate_2D = conditions[0]
        self.motion_randomness = conditions[1]
        self.learning_rate_3D = conditions[2]
        self.reconstruction_conditions = reconstruction_conditions.copy()
        self.NOS = reconstruction_conditions[2]
        self.NOS_per_section = reconstruction_conditions[4]
        self.reconstruction_conditions[2] = self.NOS_per_section + 1

        # # Pack conditions into a list
        # conditions = [noise, delta_T, NOS, theta_degrees, NOS_per_section, SRD, RDD,method,dataPiling]



    # input format, list of tuple of two elements (x,y)
    # use default again initlized in the constructor
    def append(self, snapshot):
        if self.current_snapShotIndex not in self.shotData:
            self.shotData[self.current_snapShotIndex] = snapshot
        
        if self.current_snapShotIndex not in self.shotData_original:
            self.shotData_original[self.current_snapShotIndex] = snapshot
        
        if self.current_snapShotIndex == 0:
            self.save_initial_particles(snapshot)
        else:
            self.match_previous_particle_to_current(snapshot)

        self.current_snapShotIndex += 1


# Particle ID search and assignment
    def assign_particle_id(self):
        """
        Assigns a unique ID to a particle and returns the ID.

        Returns:
        int: The unique ID assigned to the particle.
        """
        returnID = self.particle_id
        self.particle_id += 1
        return returnID

    def find_closest_particle(self, particle, shot, closest_rank=1):
            """
            Finds the index of the closest particle to a given shot.

            Args:
                particle (numpy.ndarray): An array of particle positions.
                shot (numpy.ndarray): An array of shot positions.
                closest_rank (int, optional): The rank of the closest particle to find. Defaults to 1.

            Returns:
                int: The index of the closest particle.
            """
            distances = np.linalg.norm(shot - particle, axis=1)
            minDistanceIndex = np.argmin(distances)
            if closest_rank > 1:
                for i in range(closest_rank - 1):
                    minDistanceIndex = np.argmin(np.delete(distances, minDistanceIndex))
            return minDistanceIndex
    
    def get_particle_id(self, particle, snapshotID, closest_rank=1, tolerance=0.01):
            """
            Finds the ID of the particle closest to the given particle in the specified snapshot.

            Args:
            particle (numpy.ndarray): The particle to find the closest match for.
            snapshotID (int): The ID of the snapshot to search for the particle in.
            closest_rank (int): The rank of the closest particle to return. Default is 1.
            tolerance (float): The maximum distance between the particles for them to be considered a match. Default is 0.01.

            Returns:
            tuple: The ID of the closest particle and its coordinates in the specified snapshot.

            Raises:
            KeyError: If the particle ID is not found.
            """
            target_snapshot = self.shotData[snapshotID]
            closest_particle_id_in_shot = self.find_closest_particle(particle, np.array(self.shotData[snapshotID]), closest_rank)

            for particle_id in self.particleData_2D:
                particle_relative_shotID = self.find_relative_snapshotIndex(particle_id, snapshotID)
                distance = np.linalg.norm(target_snapshot[closest_particle_id_in_shot] - self.particleData_2D[particle_id]['coords'][particle_relative_shotID])
                if distance < tolerance:
                    return particle_id, self.particleData_2D[particle_id]['coords'][particle_relative_shotID]

            print("not found")
            KeyError("particle_id not found")
    
    def get_particle_id_from_available_ids(self, particle, snapshotID, id_list):
            """
            Finds the closest particle ID in a given list of IDs to a given particle, based on their coordinates in a specific snapshot.

            Args:
                particle (int): The ID of the particle to find the closest match for.
                snapshotID (int): The index of the snapshot to use for comparing particle coordinates.
                id_list (list): A list of particle IDs to search for a match in.

            Returns:
                tuple: A tuple containing the ID of the closest matching particle and its coordinates in the specified snapshot.
            """
            target_coordinates = []
            
            for particle_id in id_list:
                particle_relative_shotID = self.find_relative_snapshotIndex(particle_id, snapshotID)
                target_coordinates.append(self.particleData_2D[particle_id]['coords'][particle_relative_shotID])
            closest_particle_id_in_shot = id_list[self.find_closest_particle(particle, np.array(target_coordinates))]
            return closest_particle_id_in_shot, self.particleData_2D[closest_particle_id_in_shot]['coords'][particle_relative_shotID]



    def get_particle_id_from_unmatched_ids(self, particle, snapshotID, matched_id_list):
            """
            Finds the closest particle ID in a given snapshot that has not been matched with any other particle IDs.

            Args:
                particle (dict): Dictionary containing the coordinates of the particle to match.
                snapshotID (int): The snapshot index to search for the closest particle.
                matched_id_list (list): List of particle IDs that have already been matched with other particles.

            Returns:
                tuple: A tuple containing the closest particle ID and its coordinates in the given snapshot.
            """
            id_list = list(range(0, len(self.particleData_2D)))

            for id in matched_id_list:
                id_list.remove(id)

            target_coordinates = []
            for particle_id in id_list:
                particle_relative_shotID = self.find_relative_snapshotIndex(particle_id, snapshotID)
                target_coordinates.append(self.particleData_2D[particle_id]['coords'][particle_relative_shotID])
            closest_particle_id_in_shot = id_list[self.find_closest_particle(particle, np.array(target_coordinates))]
            
            return closest_particle_id_in_shot, self.particleData_2D[closest_particle_id_in_shot]['coords'][particle_relative_shotID]

    
    
    def save_initial_particles(self, snapshot):
            """
            Saves the initial particles in the particleData_2D dictionary along with their coordinates and snapshot index.

            Args:
                snapshot (list): A list of particles in the current snapshot.

            Returns:
                None
            """
            for particle in snapshot:
                particle_id = self.assign_particle_id()
                if particle_id not in self.particleData_2D:
                    self.particleData_2D[particle_id] = {'coords': [], 'snapshotIndexList': [], 'snapshotIndexSet': set()}
                    
                self.particleData_2D[particle_id]['coords'].append(particle)
                self.particleData_2D[particle_id]['snapshotIndexList'].append(self.current_snapShotIndex)
                self.particleData_2D[particle_id]['snapshotIndexSet'].add(self.current_snapShotIndex)




    def rank_particle_distances(self, previous_snapshot, current_snapshot, search_radius):
                    """
                    Ranks the particles in the current snapshot based on their distance to the particles in the previous snapshot.

                    Args:
                    - previous_snapshot: A list of tuples representing the x and y coordinates of particles in the previous snapshot.
                    - current_snapshot: A list of tuples representing the x and y coordinates of particles in the current snapshot.
                    - search_radius: The maximum distance between two particles for them to be considered "close".

                    Returns:
                    - A heap queue containing tuples of the form (distance, current_particle_index, previous_particle_index), where
                        distance is the distance between the current particle and the closest particle in the previous snapshot, and
                        current_particle_index and previous_particle_index are the indices of the current particle and the closest particle
                        in the previous snapshot, respectively.
                    """
                    # turn previous snapshot into a N by 2 matrix
                    previous_snapshot = np.array(previous_snapshot).reshape(-1,2)

                    # grab the x axis
                    previous_x = previous_snapshot[:,0]
                    # grab the y axis
                    previous_y = previous_snapshot[:,1]

                    # turn current snapshot into a N by 2 matrix
                    current_snapshot = np.array(current_snapshot).reshape(-1,2)

                    # grab the x axis
                    current_x = current_snapshot[:,0]
                    # grab the y axis
                    current_y = current_snapshot[:,1]

                    # get the large matrix to fix the case where there are different number of particles in the previous and current snapshot
                    x_large_matrix = np.tile(previous_x, (len(current_x), 1)) - np.tile(current_x, (len(previous_x), 1)).T
                    y_large_matrix = np.tile(previous_y, (len(current_y), 1)) - np.tile(current_y, (len(previous_y), 1)).T

                    # print("x_large_matrix: ",x_large_matrix)
                    # print("y_large_matrix: ",y_large_matrix)

                    # get the distance matrix
                    distance_matrix = np.sqrt(x_large_matrix**2 + y_large_matrix**2)
                    # print("distance_matrix: ",distance_matrix)

                    # create a heap queue to store the ranked particles
                    ranked_particle_heapq = []

                    # Get the sorted indices of the flattened distance_matrix
                    sorted_indices = np.argsort(distance_matrix.ravel())

                    # Convert the flattened indices to 2D row and column indices
                    row_col_indices = np.unravel_index(sorted_indices, distance_matrix.shape)

                    # Print the values in distance_matrix in ascending order along with their row and column indices
                    for i in range(len(sorted_indices)):
                            row, col = row_col_indices[0][i], row_col_indices[1][i]
                            # print(f"Value: {distance_matrix[row, col]}, Row: {row}, Col: {col}")
                            distance = distance_matrix[row, col]
                            # store the closest particles together. row is the index of the curret particle, col is the index of the previous particle
                            heapq.heappush(ranked_particle_heapq, (distance, row, col))
                            # print("original ranked_particle_heapq: ",len(ranked_particle_heapq))

                    return ranked_particle_heapq
            


    def match_previous_particle_to_current(self, current_shot):
        
        """
        Matches particles from the previous snapshot to the current snapshot based on their distances.
        If the motion randomness is too high, it conducts a reconstruction to verify if the point is valid.
        If there is not enough data to do the reconstruction, it uses 2D historical velocity instead.
        
        Args:
        - current_shot: a numpy array of shape (n, 2) representing the current snapshot of particles
        
        Returns:
        - None
        """


        def is_motion_random(historical_vel, observed_vel, motion_randomness):
            """
            Determines if the motion of a particle is random based on its historical velocity and observed velocity.

            Args:
                historical_vel (float or numpy.ndarray): The historical velocity of the particle.
                observed_vel (float or numpy.ndarray): The observed velocity of the particle.
                motion_randomness (float): The threshold for determining if the motion is random.

            Returns:
                bool: True if the motion is random, False otherwise.
            """
            threshold = 1e-6
            is_random = False
            if np.isscalar(historical_vel):
                length_of_historical_vel = 1
            else:
                length_of_historical_vel = len(historical_vel)

            for i in range(length_of_historical_vel):
                if np.isscalar(observed_vel):
                    compare_observe_vel = observed_vel
                else:
                    compare_observe_vel = observed_vel[i]

                if np.isscalar(historical_vel):
                    compare_historical_vel = historical_vel
                else:
                    compare_historical_vel = historical_vel[i]

                if compare_historical_vel <= threshold:
                    # stationary, no division
                    x = abs(compare_historical_vel - compare_observe_vel) > motion_randomness
                else:
                    x = abs(compare_historical_vel - compare_observe_vel) > motion_randomness

                if x:
                    is_random = True
                    break

            return is_random


        previous_shot = self.shotData[self.current_snapShotIndex - 1]

        # create defensive copies of the previous and current shots so we can delete items to keep track without affecting the original data
        previous_shot_remain = previous_shot.copy()
        current_shot_remain = current_shot.copy()
    
        # get the ranked list of particles
        ranked_particle_list = self.rank_particle_distances(previous_shot_remain, current_shot_remain, search_radius=10)

        matched_particles_id = []


        # while there are still particles unmatched, we keep matching
        while len(previous_shot_remain) > 0 and len(current_shot_remain) > 0:
        

            # get the closest particle
            closest_particles = heapq.heappop(ranked_particle_list)
       
            previous_index = closest_particles[2]
            current_index = closest_particles[1]
            current_particle_to_match = current_shot[current_index]
            print(" ")
            print("---------------------------------")
            print("current_snapshot: ",self.current_snapShotIndex)
            # [previous_particle_id, prev_particle_coor] = self.get_particle_id(previous_shot[previous_index], self.current_snapShotIndex - 1)
            # print("particleData: ",self.particleData)
            if matched_particles_id is None:
                [previous_particle_id, prev_particle_coor] = self.get_particle_id(previous_shot[previous_index], self.current_snapShotIndex - 1)
            else:
                [previous_particle_id, prev_particle_coor] = self.get_particle_id_from_unmatched_ids(previous_shot[previous_index], self.current_snapShotIndex - 1,matched_particles_id)
            # delete the points that are matched from the defensive copies
            if not self.find_array_in_list(prev_particle_coor, previous_shot_remain) or not self.find_array_in_list(current_particle_to_match, current_shot_remain):
                # print("skip this occurance \n \n")
                continue
            
            # delete the points that are matched from the defensive copies
            for idx, particle in enumerate(previous_shot_remain):
                if np.array_equal(particle, prev_particle_coor):
                    del previous_shot_remain[idx]

                    break
            # previous_shot_remain.remove(prev_particle_coor)
          
            # explicit loop to remove the element from the current shot list
            for idx, particle in enumerate(current_shot_remain):
                if np.array_equal(particle, current_particle_to_match):
                    del current_shot_remain[idx]
                    break
            print("we are matching: ",prev_particle_coor, " with ", current_particle_to_match, " particle id: ",previous_particle_id)


            if previous_particle_id not in self.particleData_2D:
                self.add_to_particle_data(previous_particle_id)

            if self.current_snapShotIndex > 10:
                
                estimated_vel_from_historical_velocity = self.historicalLinearVelocity(previous_particle_id, 10, discount_factor=0.9)
                
                observed_vel = current_shot[current_index] - previous_shot[previous_index]
                
                if is_motion_random(estimated_vel_from_historical_velocity, observed_vel, self.motion_randomness):

                    print("motion randomness detected in particle: ",previous_particle_id)
                    # motion randomness is too high, we then conduct a reconstruction to verify if the point is valid
                    # 3D reconstruction and re sorting strategy 
                    previous_relative_index = self.find_relative_snapshotIndex(previous_particle_id, self.current_snapShotIndex - 1)

                    compensation_3D = False
                    # make sure we have enough data to do the reconstruction
                    if previous_relative_index - self.NOS_per_section > 8 and compensation_3D == True:

                        # take the previous NOS_per_section data to do the reconstruction for comparison

                        previous_2D_shots_selected = self.particleData_2D[previous_particle_id]['coords'][previous_relative_index - self.NOS_per_section : previous_relative_index]
                        print("length of previous_2D_shots_selected: ",len(previous_2D_shots_selected))

                        self.reconstruction_conditions[2] = self.NOS
                        previous_estimated_positions_single = Phase4_trace_3d(self.reconstruction_conditions, np.array(self.particleData_2D[previous_particle_id]['coords'][:]))

                        self.reconstruction_conditions[2] = self.NOS_per_section + 1
                        # take the previous NOS_per_section data and the current shot to be added to do the reconstruction and compare
                        current_estimated_positions_single = Phase4_trace_3d(self.reconstruction_conditions, np.row_stack([np.array(self.particleData_2D[previous_particle_id]['coords'][previous_relative_index - self.NOS_per_section + 1 : previous_relative_index + 1]), current_particle_to_match]))

                        learning_rate_3D = self.learning_rate_3D
                        exploitation_rate_3D = 1 - learning_rate_3D

                        # get observation velocity from the reconstruction
                        observed_vel_3D = current_estimated_positions_single[-1,:] - current_estimated_positions_single[-2,:]

                        # take the last 10 data points to do the historical velocity estimation
                        # take the last 10 rows and all columns
                        historical_vel_3D = np.mean(np.diff(previous_estimated_positions_single[-10:,:], axis=0))
        
                        # re-calculate the final position if adjustment is needed
                        adjusted_vel_3D = observed_vel_3D
                        j = 0
                        for vel in adjusted_vel_3D:
                            if is_motion_random(historical_vel_3D, observed_vel_3D[j], self.motion_randomness):
                                adjusted_vel_3D[j] = historical_vel_3D * exploitation_rate_3D + vel * learning_rate_3D

                            j += 1

                        adjusted_position_3D = current_estimated_positions_single[-1] + adjusted_vel_3D
                        adjusted_position_2D = self.particle_projection(self.alpha, adjusted_position_3D)

                        # # now update the estimated position
                        # if previous_particle_id not in self.particleData_3D:
                        #     self.particleData_3D[previous_particle_id] = {'coords': [], 'snapshotIndexList': [], 'snapshotIndexSet': set()}
                        #     self.particleData_3D[previous_particle_id]['coords'].append(current_estimated_positions_single[:-1])
                        #     self.particleData_3D[previous_particle_id]['snapshotIndexList'] = self.particleData_2D[previous_particle_id]['snapshotIndexList']
                        #     self.particleData_3D[previous_particle_id]['snapshotIndexSet'] = self.particleData_2D[previous_particle_id]['snapshotIndexSet']
                            
                            
                        # self.particleData_3D[previous_particle_id]['coords'].append(current_estimated_positions_single[-1] + adjusted_vel_3D)

                        # # add the new positions with the old, and then take average
                        # new_positions[:len(last_positions), :] = (new_positions[:len(last_positions), :] + last_positions) / 2
                        # positions_predicted[proj_used_index:proj_used_index + new_positions.shape[0], :] = new_positions
                        print("adjusted_position_2D: ",adjusted_position_2D)
                        final_xy = adjusted_position_2D

 

                    else:
                        print("not enough data to do the reconstruction, use 2D historical velocity instead")
                        # calculate position from learning factor
                        learning_rate_2D = self.learning_rate_2D
                        exploitation_rate_2D = 1 - learning_rate_2D
                        final_xy = estimated_vel_from_historical_velocity*exploitation_rate_2D + observed_vel*learning_rate_2D + previous_shot[previous_index]
                        print("final_xy: ", final_xy)
                    
                    self.particleData_2D[previous_particle_id]['coords'].append(final_xy)
                        # modify shot data to keep the consistency
                    self.shotData[self.current_snapShotIndex][current_index] = final_xy
                
                else:
                    self.particleData_2D[previous_particle_id]['coords'].append(current_particle_to_match)
                
            else:
                self.particleData_2D[previous_particle_id]['coords'].append(current_particle_to_match)

    
            self.particleData_2D[previous_particle_id]['snapshotIndexList'].append(self.current_snapShotIndex)
            
            self.particleData_2D[previous_particle_id]['snapshotIndexSet'].add(self.current_snapShotIndex)

            matched_particles_id.append(previous_particle_id)


 
        
        # print("Entered compensation mode")
        # if the new snapshot has more particles than the previous one by comparing the length of the remaining particles in the defensive copies
        if len(current_shot_remain) > len(previous_shot_remain):
            print("Entered compensation mode, more current particles than previous")
            # print("current_shot_remain: ",current_shot_remain)
            # create new unique particles and save them 
            for particle in current_shot_remain:
                previous_particle_id = self.assign_particle_id()
                if previous_particle_id not in self.particleData_2D:
                    self.particleData_2D[previous_particle_id] = {'coords': [], 'snapshotIndexList': [], 'snapshotIndexSet': set()}
                    
                self.particleData_2D[previous_particle_id]['coords'].append(particle)
                self.particleData_2D[previous_particle_id]['snapshotIndexList'].append(self.current_snapShotIndex)
                self.particleData_2D[previous_particle_id]['snapshotIndexSet'].add(self.current_snapShotIndex)

        # if the new snapshot has less particles than the previous one
        elif len(current_shot_remain) < len(previous_shot_remain):
            print("Entered compensation mode: more previous particles than current")
            # print("previous_shot_remain: ",previous_shot_remain)
            # we estimate the unmatched particle with the trajectory of the closest neighbor (current snapshot position - previous snapshot position)
            
            for prev_particle in previous_shot_remain:

                previous_particle_id, prev_particle_coor = self.get_particle_id_from_unmatched_ids(prev_particle_coor, self.current_snapShotIndex-1,matched_particles_id)
                
                # neighbor strategy
                if self.current_snapShotIndex <= 10:
                    # print("id list: ",matched_particles_id)
                    closest_neighbor_particle_id, closest_neighbor_previous_xy = self.get_particle_id_from_available_ids(prev_particle, self.current_snapShotIndex - 1, matched_particles_id)
                    # print("closest_neighbor_particle_id: ",closest_neighbor_particle_id)
                    relativeIndex = self.find_relative_snapshotIndex(closest_neighbor_particle_id, self.current_snapShotIndex)
                    # print(relativeIndex)
                    closest_neighbor_current_xy = self.particleData_2D[closest_neighbor_particle_id]['coords'][relativeIndex]
                    # print("closest_neighbor_current_xy: ",closest_neighbor_current_xy)
                    # closest_neighbor_previous_xy = self.get_coordinates_by_snapshot(closest_neighbor_particle_id, self.current_snapShotIndex - 1)
                    # Calculate the difference between current and previous coordinates (c-p)
                    difference_xy = np.array(closest_neighbor_current_xy) - np.array(closest_neighbor_previous_xy)

                    estiamted_xy = tuple(np.array(prev_particle_coor) + difference_xy)



                else:
                    previous_xy = self.particleData_2D[previous_particle_id]['coords'][self.find_relative_snapshotIndex(previous_particle_id, self.current_snapShotIndex - 1)]
                    # historical velocity strategy
                    estiamted_xy = self.historicalLinearVelocity(previous_particle_id,10,0.8) + previous_xy

                self.particleData_2D[previous_particle_id]['coords'].append(estiamted_xy)
                self.particleData_2D[previous_particle_id]['snapshotIndexList'].append(self.current_snapShotIndex)
                self.particleData_2D[previous_particle_id]['snapshotIndexSet'].add(self.current_snapShotIndex)
                self.shotData[self.current_snapShotIndex].append(np.array(estiamted_xy))

    def add_to_particle_data(self, previous_particle_id):
        self.particleData_2D[previous_particle_id] = {'coords': [], 'snapshotIndexList': [], 'snapshotIndexSet': set(), 'evaluation_state': int}


    # return historical linear velocity
    def historicalLinearVelocity(self, previous_particle_id, num_of_snapshots_to_check, discount_factor=1, direction="backward"):
        last10Coordiantes = []

        for i in range (1,num_of_snapshots_to_check + 1):
            last10Coordiantes.append(self.particleData_2D[previous_particle_id]['coords'][self.find_relative_snapshotIndex(previous_particle_id, self.current_snapShotIndex - i)])
                        # k = k**2
                        
        velocity_array = np.diff(np.array(last10Coordiantes), axis=0)
        final_velocity = np.zeros(2)
        k = 0

        for i in range(len(velocity_array)):
            if direction == "backward":
                final_velocity += velocity_array[i] * discount_factor**(i+1)
            elif direction == "forward":
                final_velocity += velocity_array[len(velocity_array)-1-i] * discount_factor**(i+1)
            else:
                # assume backward
                final_velocity += velocity_array[i] * discount_factor**(i+1)
            k += discount_factor**(i+1)
        return final_velocity/k



    def find_array_in_list(self,target, list_of_arrays):
        for idx, arr in enumerate(list_of_arrays):
            # print("equal between: ",target, " and ", arr, " is: ",np.array_equal(target, arr))
            if np.array_equal(target, arr):

                return True
        return False


    def get_coordinates_by_snapshot(self, particle_id, snapshot_index):
        if particle_id in self.particleData_2D:
            data = self.particleData_2D[particle_id]
            
            # Check if the given snapshot index exists in the list of snapshot indices
            if snapshot_index in data['snapshotIndexList']:
                # Find the index of the given snapshot index in the list
                index = data['snapshotIndexList'].index(snapshot_index)
                
                # Use the index to access the coordinates
                coordinates = data['coords'][index]
                return coordinates

        # Return None if the particle or snapshot index is not found
        return None
    
    def find_relative_snapshotIndex(self, particle_id, snapshot_index):
        if particle_id in self.particleData_2D:
            data = self.particleData_2D[particle_id]

            # print("data with particle id, ",particle_id, " is: ",data)
            # print("snapshot_index: ",snapshot_index)
            # if self.current_snapShotIndex > 141:
            #     print(self.particleData_2D)
            
            # Check if the given snapshot index exists in the list of snapshot indices
            if snapshot_index in data['snapshotIndexList']:
                # Find the index of the given snapshot index in the list
                index = data['snapshotIndexList'].index(snapshot_index)
                

                return index

        # Return None if the particle or snapshot index is not found
        KeyError("particle_id or snapshot_index not found")
        return None

    
        
    def extract_matrices(self):
        matrices = []

        for particle_id, data in self.particleData_2D.items():
            # Get the list of snapshot indices for this particle
            snapshot_indices = data['snapshotIndexList']
            
            # Initialize a matrix for this particle
            matrix = np.zeros((len(snapshot_indices), 2))  # 2 columns for x and y
            
            # Fill in the matrix with particle coordinates
            for i, snapshot_index in enumerate(snapshot_indices):
                snapshot = self.shotData[snapshot_index]
                particle_index = snapshot.index(data['coords'][i])
                x, y = snapshot[particle_index:particle_index+2]
                matrix[i] = [x, y]
            
            matrices.append(matrix)

        return matrices
    
    def get_particle_data(self):
        return self.particleData_2D
    
    def particle_projection(self, alpha, r_0):
        r_0_rotated=T(r_0,alpha)
        _, _, _, _, _, SRD, RDD,_,_ = self.reconstruction_conditions
        M_p = (SRD + RDD) / (SRD + r_0_rotated[1])
        
    
        return np.array([M_p * r_0_rotated[0], M_p * r_0_rotated[2]])
    
    # def rejectList(self, particle_id_list): 


    # def reject(self, particle_id, )
        
    # def resort()
    


    
