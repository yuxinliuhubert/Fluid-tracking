import numpy as np
from T import T  # Assuming this function rotates a particle by an angle 'alpha'
import heapq
# make sure we make enough defensive copies of the data
class ParticlePathFinder:

    def __init__(self, alpha) -> None:
        self.alpha = alpha
        self.current_snapShotIndex = 0
        self.particle_id = 0
        self.particleData = {}
        self.shotData = {}


    # input format, list of tuple of two elements (x,y)
    # use default again initlized in the constructor
    def append(self, snapshot):
        # snapshot = snapshot.tolist()
        # store the snapshot in a time sequence dictionary
        print("current_snapShotIndex: ",self.current_snapShotIndex)
        print("current shot: ",snapshot)
        if self.current_snapShotIndex != 0:
            print("previous shot: ", self.shotData[self.current_snapShotIndex - 1])

        
        if self.current_snapShotIndex not in self.shotData:
            self.shotData[self.current_snapShotIndex] = snapshot
        
        if self.current_snapShotIndex == 0:
            self.save_initial_particles(snapshot)
        else:
            # match previous particles to current

            self.match_previous_particle_to_current(snapshot)

        print("particleData: ",self.particleData)

        self.current_snapShotIndex += 1

   


    # def appendCustomAngle(self, snapshot, alpha):
    #     if self.current_snapShotIndex == 0:
    #         self.outputMatrix = snapshot
    #         return
    #     self.adjust_matrix_size(snapshot, alpha)
    #     self.current_snapShotIndex += 1
    
        # newShot
        # self.outputMatrix = np.concatenate((self.outputMatrix, self.newShot), axis=1)
        # get the rows and columns of the data set
        # compare column of data set with existing matrix, if N_dataset > matrix, expand the matrix, otherwise:
        # 1. match the existing points with current data set 
        # 2. for the remaining points in the existing matrix, "hold the position" until more data is avaialble 

        # data set shape -> 1 row, N*2 columns, where N is the number of particles. 1st is x1, 2nd column is y1, 3rd is x2, 4th is y2, etc.
    
    # def adjust_matrix_size(self, snapshot, alpha):
    #     n_particles = snapshot.shape[1] // 2
    #      # If the number of particles in the dataset is greater than in the output matrix,
    #     # then expand the matrix
    #     if snapshot.shape[1] > self.outputMatrix.shape[1]:
    #         num_particles_to_match = self.outputMatrix.shape[1]
    #         additional_columns = 2 * n_particles - self.outputMatrix.shape[1]
    #         self.outputMatrix = np.concatenate((self.outputMatrix, np.zeros((self.outputMatrix.shape[0], additional_columns))), axis=1)            
    #         self.match_particles(snapshot, alpha, num_particles_to_match)

    #     else:
    #         self.match_particles(snapshot, alpha, self.outputMatrix.shape[1] // 2)

        


# find the closest neighbors of the unmatched particles with matched particles in the same shot. translate the movement of the matched particle to unmatched
# 2. in cases where the unmatched parcile has a history, we can just use the historical movement to predict the next position.
    # def match_particles(self, current_snapshot, alpha, num_particles_to_match):
    #     if self.current_snapShotIndex == 0:
    #         print("You cannot track particles without at least 2 snapshots.")
    #         return
        
    #     # get the previous snapshot
    #     previous_snapshot = self.outputMatrix[self.current_snapShotIndex, :]
    #     prev_N = previous_snapshot.shape[1] // 2
    #     # N = number of particles in the current snapshot
    #     M, N = current_snapshot.shape[0], current_snapshot.shape[1] // 2  # Assuming each particle has 2 coordinates (x, z)
    #     matched_particles = np.zeros((1, N * 2))

    #     # match the available particles (previous and current both available)
    #     for n in range(0, 2 * num_particles_to_match, 2):

    #         # if there are no more particles in the current snapshot, we estimate the positions of the unmatched particles using neighboring particles'
    #         # movements from previous to current snapshot
    #         if current_snapshot.shape[1] == 0:

    #             continue
    #         particle = previous_snapshot[n:n+2]
    #         closest_particle_index = self.find_closest_particle(particle, current_snapshot)
    #         matched_particles[2*closest_particle_index:2*closest_particle_index+2] = particle
    #         # remove the matched particle from the previous snapshot
    #         current_snapshot = np.delete(current_snapshot, closest_particle_index, axis=1)

    #     # match the remaining particles (only current available)
    #     if current_snapshot.shape[1] > 1:
    #         matched_particles = np.concatenate((matched_particles, current_snapshot), axis=1)



    #     self.outputMatrix = np.concatenate((self.outputMatrix, matched_particles), axis=0)
    #     return 
    

    # def estimate_position_with_neighbor(self, particle, shot):
    #     # find the closest neighbor of the particle in the current shot, and use the movement of the neighbor to predict the position of the particle



    def find_closest_particle(self, particle, shot, closest_rank=1):
        print("shot in find closest particle: ",shot)
        print("particle in find closest particle: ",particle)
        # print("shotindex: ",self.current_snapShotIndex)
        distances = np.linalg.norm(shot - particle, axis=1)
        print("distances: ", distances)
        print("the argmin is: ",np.argmin(distances))

        minDistanceIndex = np.argmin(distances)

        if closest_rank > 1:
            for i in range(closest_rank - 1):
                minDistanceIndex = np.argmin(np.delete(distances, minDistanceIndex))
        return minDistanceIndex
    
    def assign_particle_id(self):
        returnID = self.particle_id
        self.particle_id += 1
        return returnID

    def get_particle_id(self, particle,snapshotID, closest_rank=1):
        target_snapshot = self.shotData[snapshotID]
        closest_particle_id_in_shot = self.find_closest_particle(particle, np.array(self.shotData[snapshotID]),closest_rank)

        print("particle in get particle id: ",particle)
        print("closest_particle_coor: ",closest_particle_id_in_shot)
        # print("target_snapshot: ",target_snapshot)
        print("particleData with id: ",self.particleData[closest_particle_id_in_shot])
        # print("snapshotID relative: ",particle_relative_shotID)
        # print("particleData on this shot: ",self.particleData[closest_particle_id_in_shot]['coords'][particle_relative_shotID])
        for particle_id in self.particleData:

            print("iterating at particle_id: ",particle_id)
            # print("particleData_individual: ", self.particleData[particle_id]['coords'][snapshotID])
            # print("target_snapshot[closest_particle_id_in_shot]: ",target_snapshot[closest_particle_id_in_shot])
            # print("self.particleData[particle_id]['coords'][snapshotID]: ",self.particleData[particle_id]['coords'][snapshotID])
            particle_relative_shotID = self.find_relative_snapshotIndex(particle_id, snapshotID)
            if np.array_equal(target_snapshot[closest_particle_id_in_shot], self.particleData[particle_id]['coords'][particle_relative_shotID]):
                print("found the particle id: ",particle_id)
                # particle_relative_shotID = self.find_relative_snapshotIndex(particle_id, snapshotID)
                print("particle_relative_shotID: ",particle_relative_shotID)
                print(self.particleData[particle_id]['coords'][particle_relative_shotID])
            
                return particle_id, self.particleData[particle_id]['coords'][particle_relative_shotID]
            
        print("not found")
        KeyError("particle_id not found")

        

        # if matched we have to eliminate the ones with matching.
        # print("particleData: ", self.particleData)
        # for particle_id in self.particleData:

        #     print("particle_id: ",particle_id)
        #     print("particleData_individual: ", self.particleData[particle_id]['coords'][snapshotID])
        #     # print("target_snapshot[closest_particle_coor]: ",target_snapshot[closest_particle_coor])
        #     if target_snapshot[closest_particle_coor].all() == self.particleData[particle_id]['coords'][snapshotID].all():

        #         return particle_id,self.particleData[particle_id]['coords'][snapshotID]
    
    def save_initial_particles(self, snapshot):
        for particle in snapshot:
            particle_id = self.assign_particle_id()
            # if self.paricleData is None:
            #     self.paricleData = {particle_id: particle}
            if particle_id not in self.particleData:
                self.particleData[particle_id] = {'coords': [], 'snapshotIndexList': [], 'snapshotIndexSet': set()}
                
            self.particleData[particle_id]['coords'].append(particle)
            self.particleData[particle_id]['snapshotIndexList'].append(self.current_snapShotIndex)
            self.particleData[particle_id]['snapshotIndexSet'].add(self.current_snapShotIndex)




    def rank_particle_distances(self, previous_snapshot, current_snapshot, search_radius):
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
        previous_shot = self.shotData[self.current_snapShotIndex - 1]

        # create defensive copies of the previous and current shots so we can delete items to keep track without affecting the original data
        previous_shot_remain = previous_shot.copy()
        
        current_shot_remain = current_shot.copy()

        

        # get the ranked list of particles
        ranked_particle_list = self.rank_particle_distances(previous_shot_remain, current_shot_remain, search_radius=10)



        # while there are still particles unmatched, we keep matching
        while len(previous_shot_remain) > 0 and len(current_shot_remain) > 0:
            

            # get the closest current and previous particle in the ranked list

            
            # particle = current_shot_remain[0]
            # print("previous shot list: ",previous_shot_remain)
            # print("current shot list: ",current_shot_remain)

    

            # get the closest particle
            closest_particles = heapq.heappop(ranked_particle_list)
            # print("ranked_particle_heapq: ",ranked_particle_list)
            # print("size of ranked_particle_list:" ,len(ranked_particle_list))
            # print("size of previous_shot_remain:" ,len(previous_shot_remain))
            # print("size of current_shot_remain:" ,len(current_shot_remain))
            # print("closest_particles: ",closest_particles)
            previous_index = closest_particles[2]
            current_index = closest_particles[1]
            current_particle_to_match = current_shot[current_index]

            
            # print("selected current: ",current_particle_to_match)

            [previous_particle_id, prev_particle_coor] = self.get_particle_id(previous_shot[previous_index], self.current_snapShotIndex - 1)
            print("selected prev: ",previous_shot[previous_index])
            print("previous_particle_id put into particleData: ",self.particleData[previous_particle_id])
            print("particleData: ",self.particleData)
            
            # delete the points that are matched from the defensive copies
            tuple_to_find = tuple(prev_particle_coor)
            
            # print("compare previous_shot_coor: ", previous_shot[previous_index], " with prev_particle_coor",  prev_particle_coor)
            # print("previous_shot_remain element type,", type(previous_shot_remain[0]))
            # print("prev coor found: ",self.find_array_in_list(prev_particle_coor, previous_shot_remain))
            # print("current coor found: ",self.find_array_in_list(current_particle_to_match, current_shot_remain))
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
            print("we are matching: ",prev_particle_coor, " with ", current_particle_to_match)

            if previous_particle_id not in self.particleData:
                self.particleData[previous_particle_id] = {'coords': [], 'snapshotIndexList': [], 'snapshotIndexSet': set()}
                
            self.particleData[previous_particle_id]['coords'].append(current_particle_to_match)
            self.particleData[previous_particle_id]['snapshotIndexList'].append(self.current_snapShotIndex)
            self.particleData[previous_particle_id]['snapshotIndexSet'].add(self.current_snapShotIndex)

            print("------------------")

        print("Entered compensation mode")
        # if the new snapshot has more particles than the previous one by comparing the length of the remaining particles in the defensive copies
        if len(current_shot_remain) > len(previous_shot_remain):
            print("current_shot_remain: ",current_shot_remain)
            # create new unique particles and save them 
            for particle in current_shot_remain:
                previous_particle_id = self.assign_particle_id()
                if previous_particle_id not in self.particleData:
                    self.particleData[previous_particle_id] = {'coords': [], 'snapshotIndexList': [], 'snapshotIndexSet': set()}
                    
                self.particleData[previous_particle_id]['coords'].append(particle)
                self.particleData[previous_particle_id]['snapshotIndexList'].append(self.current_snapShotIndex)
                self.particleData[previous_particle_id]['snapshotIndexSet'].add(self.current_snapShotIndex)

        # if the new snapshot has less particles than the previous one
        elif len(current_shot_remain) < len(previous_shot_remain):
            print("previous_shot_remain: ",previous_shot_remain)
            # we estimate the unmatched particle with the trajectory of the closest neighbor (current snapshot position - previous snapshot position) 
            for prev_particle in previous_shot_remain:

                previous_particle_id, prev_particle_coor = self.get_particle_id(prev_particle, self.current_snapShotIndex - 1)

                closest_neighbor_particle_id, closest_neighbor_previous_xy = self.get_particle_id(prev_particle_coor, self.current_snapShotIndex-1, closest_rank=2)
                # print("closest_neighbor_particle_id: ",closest_neighbor_particle_id)
                closest_neighbor_current_xy = self.particleData[closest_neighbor_particle_id]['coords'][self.current_snapShotIndex]
                # print("closest_neighbor_current_xy: ",closest_neighbor_current_xy)
                # closest_neighbor_previous_xy = self.get_coordinates_by_snapshot(closest_neighbor_particle_id, self.current_snapShotIndex - 1)
                # Calculate the difference between current and previous coordinates (c-p)
                difference_xy = np.array(closest_neighbor_current_xy) - np.array(closest_neighbor_previous_xy)
                print("difference_xy: ",difference_xy)
                estiamted_xy = tuple(np.array(prev_particle_coor) + difference_xy)
                self.particleData[previous_particle_id]['coords'].append(estiamted_xy)
                self.particleData[previous_particle_id]['snapshotIndexList'].append(self.current_snapShotIndex)
                self.particleData[previous_particle_id]['snapshotIndexSet'].add(self.current_snapShotIndex)
                self.shotData[self.current_snapShotIndex].append(np.array(estiamted_xy))


    def find_array_in_list(self,target, list_of_arrays):
        for idx, arr in enumerate(list_of_arrays):
            print("equal between: ",target, " and ", arr, " is: ",np.array_equal(target, arr))
            if np.array_equal(target, arr):

                return True
        return False


    # def get_target_snapshot(self, target_snapshotID):
    #     if target_snapshotID in self.particleData:
    #         return self.particleData[target_snapshotID]
    #     else:
    #         return None
    
    def get_coordinates_by_snapshot(self, particle_id, snapshot_index):
        if particle_id in self.particleData:
            data = self.particleData[particle_id]
            
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
        if particle_id in self.particleData:
            data = self.particleData[particle_id]
            
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

        for particle_id, data in self.particleData.items():
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
        return self.particleData


# use the existing particles that have a match to predict the ones that do not have a match. 



# class OriginalPathKeeper:
    
    
#     def __init__(self, alpha) -> None:
#         self.alpha = alpha
#         self.current_snapShotIndex = 0
#         self.particle_id = 0
#         self.paricleData = {}
#         self.shotData = {}
#         self.columnIndex = 0

#     # use default again initlized in the constructor
#     def append(self, snapshot, columnIndex):
#         self.columnIndex = columnIndex
#         # snapshot = snapshot.tolist()
#         # store the snapshot in a time sequence dictionary
#         if self.current_snapShotIndex not in self.shotData:
#             self.shotData[self.current_snapShotIndex] = snapshot
        
#         if self.current_snapShotIndex == 0:
#             self.save_initial_particles(snapshot)
#         else:
#             # match previous particles to current

#             self.match_previous_particle_to_current(snapshot)

#         self.current_snapShotIndex += 1

#     def save_initial_particles(self, snapshot):
#         for particle in snapshot:
#             particle_id = self.columnIndex
#             if particle_id not in self.particleData:
#                 self.particleData[particle_id] = {'coords': [], 'snapshotIndexList': [], 'snapshotIndexSet': set()}
                
#             self.particleData[particle_id]['coords'].append(particle)
#             self.particleData[particle_id]['snapshotIndexList'].append(self.current_snapShotIndex)
#             self.particleData[particle_id]['snapshotIndexSet'].add(self.current_snapShotIndex)

#             self.columnIndex += 1

#     def match_previous_particle_to_current(self, snapshot):
#         previous_shot = self.shotData[self.snapShotIndex - 1]

#         # create defensive copies of the previous and current shots so we can delete items to keep track without affecting the original data
#         previous_shot_remain = previous_shot.copy()
#         current_shot_remain = snapshot.copy()
#         for i in range(0, min(len(previous_shot),len(snapshot)) // 2):

#             particle = tuple(snapshot[i])
#             [particle_id, prev_particle_coor] = self.get_particle_id(particle, self.current_snapShotIndex - 1)


#             # if prev_particle_coor not in previous_shot_remain or :
#             #     continue
#             # delete the points that are matched from the defensive copies
#             previous_shot_remain.remove(prev_particle_coor)
#             current_shot_remain.remove(particle)

#             if particle_id not in self.particleData:
#                 self.particleData[particle_id] = {'coords': [], 'snapshotIndexList': [], 'snapshotIndexSet': set()}
                
#             self.particleData[particle_id]['coords'].append(particle)
#             self.particleData[particle_id]['snapshotIndexList'].append(self.snapShotIndex)
#             self.particleData[particle_id]['snapshotIndexSet'].add(self.snapShotIndex)

#         # if the new snapshot has more particles than the previous one by comparing the length of the remaining particles in the defensive copies
#         if len(current_shot_remain) > len(previous_shot_remain):
#             # create new unique particles and save them 
#             for particle in current_shot_remain:
#                 particle_id = self.assign_particle_id()
#                 if particle_id not in self.particleData:
#                     self.particleData[particle_id] = {'coords': [], 'snapshotIndexList': [], 'snapshotIndexSet': set()}
                    
#                 self.particleData[particle_id]['coords'].append(particle)
#                 self.particleData[particle_id]['snapshotIndexList'].append(self.current_snapShotIndex)
#                 self.particleData[particle_id]['snapshotIndexSet'].add(self.current_snapShotIndex)

#         # if the new snapshot has less particles than the previous one
#         elif len(current_shot_remain) < len(previous_shot_remain):
#             # we estimate the unmatched particle with the trajectory of the closest neighbor (current snapshot position - previous snapshot position) 
#             for prev_particle in previous_shot_remain:

#                 closest_neighbor_particle_id = self.get_particle_id(prev_particle_coor, self.snapShotIndex-1)
#                 closest_neighbor_current_xy = self.get_coordinates_by_snapshot(closest_neighbor_particle_id, self.snapShotIndex)
#                 closest_neighbor_previous_xy = self.get_coordinates_by_snapshot(closest_neighbor_particle_id, self.snapShotIndex - 1)
#                 # Calculate the difference between current and previous coordinates (c-p)
#                 difference_xy = np.array(closest_neighbor_current_xy) - np.array(closest_neighbor_previous_xy)
#                 estiamted_xy = tuple(np.array(particle) + difference_xy)
#                 self.particleData[particle_id]['coords'].append(estiamted_xy)
#                 self.particleData[particle_id]['snapshotIndexList'].append(self.snapShotIndex)
#                 self.particleData[particle_id]['snapshotIndexSet'].add(self.snapShotIndex)

#     def get_particle_id(self, particle,snapshotID):
#         target_snapshot = self.shotData[snapshotID]
#         closest_particle_coor = self.find_closest_particle(particle, np.array(self.shotData[snapshotID]))

#         for particle_id in self.particleData:
#             if target_snapshot[closest_particle_coor] in self.particleData[particle_id]['coords']:
#                 return tuple(particle_id,self.particleData[particle_id]['coords'])

#     def assign_particle_id(self):
#         returnID = self.particle_id
#         self.particle_id += 1
#         return returnID