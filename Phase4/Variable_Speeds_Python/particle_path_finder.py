import numpy as np
from T import T  # Assuming this function rotates a particle by an angle 'alpha'

# make sure we make enough defensive copies of the data
class ParticlePathFinder:

    def __init__(self, alpha) -> None:
        self.alpha = alpha
        self.current_snapShotIndex = 0
        self.particle_id = 0
        self.paricleData = {}
        self.shotData = {}


    # input format, list of tuple of two elements (x,y)
    # use default again initlized in the constructor
    def append(self, snapshot):
        # snapshot = snapshot.tolist()
        # store the snapshot in a time sequence dictionary
        if self.current_snapShotIndex not in self.shotData:
            self.shotData[self.current_snapShotIndex] = snapshot
        
        if self.current_snapShotIndex == 0:
            self.save_initial_particles(snapshot)
        else:
            # match previous particles to current

            self.match_previous_particle_to_current(snapshot)

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



    def find_closest_particle(self, particle, shot):
        distances = np.linalg.norm(shot - particle, axis=1)
        return np.argmin(distances)
    
    def assign_particle_id(self):
        returnID = self.particle_id
        self.particle_id += 1
        return returnID

    def get_particle_id(self, particle,snapshotID):
        target_snapshot = self.shotData[snapshotID]
        closest_particle_coor = self.find_closest_particle(particle, np.array(self.shotData[snapshotID]))

        for particle_id in self.particleData:
            if target_snapshot[closest_particle_coor] in self.particleData[particle_id]['coords']:
                return tuple(particle_id,self.particleData[particle_id]['coords'])
    
    def save_initial_particles(self, snapshot):
        for particle in snapshot:
            particle_id = self.assign_particle_id()
            if particle_id not in self.particleData:
                self.particleData[particle_id] = {'coords': [], 'snapshotIndexList': [], 'snapshotIndexSet': set()}
                
            self.particleData[particle_id]['coords'].append(particle)
            self.particleData[particle_id]['snapshotIndexList'].append(self.snapShotIndex)
            self.particleData[particle_id]['snapshotIndexSet'].add(self.snapShotIndex)

    def match_previous_particle_to_current(self, snapshot):
        previous_shot = self.shotData[self.snapShotIndex - 1]

        # create defensive copies of the previous and current shots so we can delete items to keep track without affecting the original data
        previous_shot_remain = previous_shot.copy()
        current_shot_remain = snapshot.copy()
        for i in range(0, min(len(previous_shot),len(snapshot))):

            particle = tuple(snapshot[i])
            [particle_id, prev_particle_coor] = self.get_particle_id(particle, self.snapShotIndex - 1)

            # delete the points that are matched from the defensive copies
            previous_shot_remain.remove(prev_particle_coor)
            current_shot_remain.remove(particle)

            if particle_id not in self.particleData:
                self.particleData[particle_id] = {'coords': [], 'snapshotIndexList': [], 'snapshotIndexSet': set()}
                
            self.particleData[particle_id]['coords'].append(particle)
            self.particleData[particle_id]['snapshotIndexList'].append(self.snapShotIndex)
            self.particleData[particle_id]['snapshotIndexSet'].add(self.snapShotIndex)

        # if the new snapshot has more particles than the previous one by comparing the length of the remaining particles in the defensive copies
        if len(current_shot_remain) > len(previous_shot_remain):
            # create new unique particles and save them 
            for particle in current_shot_remain:
                particle_id = self.assign_particle_id()
                if particle_id not in self.particleData:
                    self.particleData[particle_id] = {'coords': [], 'snapshotIndexList': [], 'snapshotIndexSet': set()}
                    
                self.particleData[particle_id]['coords'].append(particle)
                self.particleData[particle_id]['snapshotIndexList'].append(self.snapShotIndex)
                self.particleData[particle_id]['snapshotIndexSet'].add(self.snapShotIndex)

        # if the new snapshot has less particles than the previous one
        elif len(current_shot_remain) < len(previous_shot_remain):
            # we estimate the unmatched particle with the trajectory of the closest neighbor (current snapshot position - previous snapshot position) 
            for prev_particle in previous_shot_remain:

                closest_neighbor_particle_id = self.get_particle_id(prev_particle_coor, self.snapShotIndex-1)
                closest_neighbor_current_xy = self.get_coordinates_by_snapshot(closest_neighbor_particle_id, self.snapShotIndex)
                closest_neighbor_previous_xy = self.get_coordinates_by_snapshot(closest_neighbor_particle_id, self.snapShotIndex - 1)
                # Calculate the difference between current and previous coordinates (c-p)
                difference_xy = np.array(closest_neighbor_current_xy) - np.array(closest_neighbor_previous_xy)
                estiamted_xy = tuple(np.array(particle) + difference_xy)
                self.particleData[particle_id]['coords'].append(estiamted_xy)
                self.particleData[particle_id]['snapshotIndexList'].append(self.snapShotIndex)
                self.particleData[particle_id]['snapshotIndexSet'].add(self.snapShotIndex)


    def get_target_snapshot(self, target_snapshotID):
        if target_snapshotID in self.particleData:
            return self.particleData[target_snapshotID]
        else:
            return None
    
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



class OriginalPathKeeper:
    
    
    def __init__(self, alpha) -> None:
        self.alpha = alpha
        self.current_snapShotIndex = 0
        self.particle_id = 0
        self.paricleData = {}
        self.shotData = {}
        self.columnIndex = 0

    # use default again initlized in the constructor
    def append(self, snapshot, columnIndex):
        self.columnIndex = columnIndex
        # snapshot = snapshot.tolist()
        # store the snapshot in a time sequence dictionary
        if self.current_snapShotIndex not in self.shotData:
            self.shotData[self.current_snapShotIndex] = snapshot
        
        if self.current_snapShotIndex == 0:
            self.save_initial_particles(snapshot)
        else:
            # match previous particles to current

            self.match_previous_particle_to_current(snapshot)

        self.current_snapShotIndex += 1

    def save_initial_particles(self, snapshot):
        for particle in snapshot:
            particle_id = self.columnIndex
            if particle_id not in self.particleData:
                self.particleData[particle_id] = {'coords': [], 'snapshotIndexList': [], 'snapshotIndexSet': set()}
                
            self.particleData[particle_id]['coords'].append(particle)
            self.particleData[particle_id]['snapshotIndexList'].append(self.snapShotIndex)
            self.particleData[particle_id]['snapshotIndexSet'].add(self.snapShotIndex)

            self.columnIndex += 1

    def match_previous_particle_to_current(self, snapshot):
        previous_shot = self.shotData[self.snapShotIndex - 1]

        # create defensive copies of the previous and current shots so we can delete items to keep track without affecting the original data
        previous_shot_remain = previous_shot.copy()
        current_shot_remain = snapshot.copy()
        for i in range(0, min(len(previous_shot),len(snapshot)) // 2):

            particle = tuple(snapshot[i])
            [particle_id, prev_particle_coor] = self.get_particle_id(particle, self.snapShotIndex - 1)

            # delete the points that are matched from the defensive copies
            previous_shot_remain.remove(prev_particle_coor)
            current_shot_remain.remove(particle)

            if particle_id not in self.particleData:
                self.particleData[particle_id] = {'coords': [], 'snapshotIndexList': [], 'snapshotIndexSet': set()}
                
            self.particleData[particle_id]['coords'].append(particle)
            self.particleData[particle_id]['snapshotIndexList'].append(self.snapShotIndex)
            self.particleData[particle_id]['snapshotIndexSet'].add(self.snapShotIndex)

        # if the new snapshot has more particles than the previous one by comparing the length of the remaining particles in the defensive copies
        if len(current_shot_remain) > len(previous_shot_remain):
            # create new unique particles and save them 
            for particle in current_shot_remain:
                particle_id = self.assign_particle_id()
                if particle_id not in self.particleData:
                    self.particleData[particle_id] = {'coords': [], 'snapshotIndexList': [], 'snapshotIndexSet': set()}
                    
                self.particleData[particle_id]['coords'].append(particle)
                self.particleData[particle_id]['snapshotIndexList'].append(self.snapShotIndex)
                self.particleData[particle_id]['snapshotIndexSet'].add(self.snapShotIndex)

        # if the new snapshot has less particles than the previous one
        elif len(current_shot_remain) < len(previous_shot_remain):
            # we estimate the unmatched particle with the trajectory of the closest neighbor (current snapshot position - previous snapshot position) 
            for prev_particle in previous_shot_remain:

                closest_neighbor_particle_id = self.get_particle_id(prev_particle_coor, self.snapShotIndex-1)
                closest_neighbor_current_xy = self.get_coordinates_by_snapshot(closest_neighbor_particle_id, self.snapShotIndex)
                closest_neighbor_previous_xy = self.get_coordinates_by_snapshot(closest_neighbor_particle_id, self.snapShotIndex - 1)
                # Calculate the difference between current and previous coordinates (c-p)
                difference_xy = np.array(closest_neighbor_current_xy) - np.array(closest_neighbor_previous_xy)
                estiamted_xy = tuple(np.array(particle) + difference_xy)
                self.particleData[particle_id]['coords'].append(estiamted_xy)
                self.particleData[particle_id]['snapshotIndexList'].append(self.snapShotIndex)
                self.particleData[particle_id]['snapshotIndexSet'].add(self.snapShotIndex)

    def get_particle_id(self, particle,snapshotID):
        target_snapshot = self.shotData[snapshotID]
        closest_particle_coor = self.find_closest_particle(particle, np.array(self.shotData[snapshotID]))

        for particle_id in self.particleData:
            if target_snapshot[closest_particle_coor] in self.particleData[particle_id]['coords']:
                return tuple(particle_id,self.particleData[particle_id]['coords'])

    def assign_particle_id(self):
        returnID = self.particle_id
        self.particle_id += 1
        return returnID