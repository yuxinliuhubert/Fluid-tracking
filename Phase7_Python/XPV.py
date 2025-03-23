#!/usr/bin/env python
# coding: utf-8

# Enter parameters here

# In[ ]:


# install all the necessary packages

import subprocess
import sys

try:
    from importlib.metadata import distribution
except ImportError:
    # Fallback for Python<3.8
    from importlib_metadata import distribution

def is_package_installed(package):
    try:
        distribution(package)
        return True
    except:
        return False

def install_packages(packages, python_executable):
    for package in packages:
        if not is_package_installed(package):
            subprocess.check_call([python_executable, "-m", "pip", "install", package])
        else:
            print(f"Package '{package}' is already installed.")

if __name__ == "__main__":
    packages_to_install = ['numpy', 'pandas', 'scipy', 'plotly','sympy','openpyxl']  # Add your packages here
    python_executable = sys.executable  # Uses the Python interpreter that runs the script
    install_packages(packages_to_install, python_executable)



# numpy
import numpy as np

# scipy
from scipy import integrate
from scipy.signal import savgol_filter
from scipy.interpolate import UnivariateSpline
from scipy.stats import norm

import pandas as pd

import time
import heapq
import os
from generateTestPositions import generateTestPositions
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib
import plotly.graph_objects as go
# matplotlib.use('TkAgg')  # or another interactive backend
import random

from sympy import symbols, Eq, solve


# In[ ]:


# dataFolder = r"5"
dataFolder = r"Synthetic_1particle_linearx"
folderName = r"C:\Users\10032\Documents\GitHub\Fluid-tracking\Phase7_python\Unsorted" 
file_prefix = "NOS"

# Testing with dummy data...only required when data is not available for code testing...generates projections based on magnification, rotation and known velocity profile

NumOfDataPoints = 1
clusterness = 0.1 # smaller number the more clustered


# Input conditions
noise = 1e-3
theta_degrees = 1
# rev = 2  # revolutions of camera for the entire process
# NOS = int(rev * 360 / theta_degrees)
NOS = 360
NOS_per_section = 140  # must be larger than 5 to satisfy equations
delta_T = 0.1
# delta_T = 1/68
camera_speed = 0.5  # in Hz or revolution per second
SOD = 38  # mm, Source-Reference Distance
ODD = 462  # mm, Reference-Detector (screen) Distance
# SOD = 15  # mm, Source-Reference Distance
# ODD = 400  # mm, Reference-Detector (screen) Distance


# offset = [243.5, 97.5]
offset = [0.0,0.0]
# pixelResolution = 0.172  # every pixel is equal to mm
pixelResolution = 1  # every pixel is equal to mm
method = 'acceleration'
dataPiling = 'overlap'

# Auto-calculations of the rest of the parameters derived from the setting above
# delta_T = camera_speed * theta_degrees / 360
camera_speed = 360* delta_T / theta_degrees
shots_per_second = 1 / delta_T

# Define the velocity function
# v = lambda t: [0.9 * np.sin(t), 0.9 * np.cos(t), 1]

# AI conditions
learning_rate_2D =1
motion_randomness = 3
learning_rate_3D =0.3


# Pack conditions into a list
conditions = [noise, delta_T, NOS, theta_degrees, NOS_per_section, SOD, ODD,method,dataPiling]



# In[ ]:


# Generate test positions
# NumberOfTestPoints = 1
# initial_positions = np.zeros((NumOfDataPoints,3))
# initial_positions[0] = [6,4,8]
# v = []

# # v.append(lambda t: [3*np.sin(t), 2*np.cos(t), np.sin(t)])
# v.append(lambda t: [0.1,0.05,-0.1])




real_positions_folder = r"C:\Users\10032\Documents\GitHub\Fluid-tracking\Phase7_python\Real_positions"

print("list of files: ",os.listdir(real_positions_folder))
sorted_filenames = sorted(os.listdir(real_positions_folder), key=lambda x: int(x.split(file_prefix)[1].split('.csv')[0]))
print(len(sorted_filenames))
k = 0

for fileIndex in range(min(len(sorted_filenames), NOS)):
    file = sorted_filenames[fileIndex]
    if file.endswith(".csv"):
        filename = os.path.join(real_positions_folder, file)
        input_data = pd.read_csv(filename, header=None)
        input_data = np.array(np.transpose(input_data))
        if fileIndex == 0:
            real_positions = input_data
        else:
            real_positions = np.row_stack((real_positions, input_data))
    # print("fileIndex: ",fileIndex)
# print("real_positions: ",real_positions)
# _, columnNum = real_positions.shape
# # Create the figure and axes
# fig3 = plt.figure()
# ax3 = fig3.add_subplot(111, projection='3d')
# for i in range(int(columnNum/3)):
#     plotting_single(real_positions[:,i*3:i*3+3],i,ax3)

# plt.show()

sorted_positions_folder = r"C:\Users\10032\Documents\GitHub\Fluid-tracking\Phase7_python\Sorted"

print("list of files: ",os.listdir(sorted_positions_folder))
sorted_filenames = sorted(os.listdir(sorted_positions_folder), key=lambda x: int(x.split(file_prefix)[1].split('.csv')[0]))
print(len(sorted_filenames))
k = 0

for fileIndex in range(min(len(sorted_filenames), NOS)):
    file = sorted_filenames[fileIndex]
    if file.endswith(".csv"):
        filename = os.path.join(sorted_positions_folder, file)
        input_data = pd.read_csv(filename, header=None)
        input_data = np.array(np.transpose(input_data))
        # print("input_data: ",input_data)
        if fileIndex == 0:
            xz_proj = input_data
        else:
            xz_proj = np.row_stack((xz_proj, input_data))

NumOfDataPoints = int(len(xz_proj[0])/2)

# calculate synthetic data 
# xz_proj = np.zeros((NOS, NumOfDataPoints*2))
# real_positions = np.zeros((NOS, NumOfDataPoints*3))
# # Generate test positions
# for i in range(NumOfDataPoints):
#     # vel = v[i]
#     # xz_proj[:,i*2:i*2+2], real_positions[:,i*3:i*3+3]= generateTestPositions(vel, initial_positions[i], conditions)
#     xz_proj[:,i*2:i*2+2], _= generateTestPositions(vel, initial_positions[i], conditions)


# In[ ]:


# shared functions
def map_range(x, x_min, x_max, y_min, y_max):
    # Apply linear mapping
    scaled = y_min + (y_max - y_min) * ((x - x_min) / (x_max - x_min))

    # Apply floor and ceiling conditions
    scaled = np.where(x <= x_min, y_min, scaled)
    scaled = np.where(x >= x_max, y_max, scaled)

    return scaled

def rotation(r1, alpha):#why -alpha?


    rotation_matrix = np.array([
        [np.cos(-alpha), -np.sin(-alpha), 0],
        [np.sin(-alpha),  np.cos(-alpha), 0],
        [0,               0,              1]
    ])
    r2 = np.matmul(rotation_matrix, r1)#does order matter?
    return r2




# In[ ]:


# 2D sorting functions
class pf:
    def __init__(self, alpha,conditions, reconstruction_conditions) -> None:
        self.alpha = alpha
        self.current_snapShotIndex = 0
        self.particle_id = 0
        self.particleData_2D = {}
        self.shotData = {}
        self.original_shotData = {}
        self.learning_rate_2D = conditions[0]
        self.corrected_learning_rates_data = {}
        self.motion_randomness = conditions[1]
        self.learning_rate_3D = conditions[2]
        self.reconstruction_conditions = reconstruction_conditions.copy()
        self.NOS = reconstruction_conditions[2]
        self.NOS_per_section = reconstruction_conditions[4]
        self.reconstruction_conditions[2] = self.NOS_per_section + 1
        self.learning_rate_corrected = False
        self.SOD = reconstruction_conditions[5]
        self.RDD = reconstruction_conditions[6]

        # # Pack conditions into a list
        # conditions = [noise, delta_T, NOS, theta_degrees, NOS_per_section, SRD, RDD,method,dataPiling]
     # input format, list of tuple of two elements (x,y)

    def assign_particle_id(self):
        returnID = self.particle_id
        self.particle_id += 1
        return returnID
    
    def append(self, snapshot):
        # snapshot = snapshot.tolist()
        # store the snapshot in a time sequence dictionary
        # print("current_snapShotIndex: ",self.current_snapShotIndex)
        # print("current shot: ",snapshot)
        # if self.current_snapShotIndex != 0:
        #     print("previous shot: ", self.shotData[self.current_snapShotIndex - 1])

    
        if self.current_snapShotIndex not in self.shotData:
            self.shotData[self.current_snapShotIndex] = snapshot

        if self.current_snapShotIndex not in self.original_shotData:
            self.original_shotData[self.current_snapShotIndex] = snapshot.copy()
        
        
        if self.current_snapShotIndex == 0:
            self.save_initial_particles(snapshot)
        else:
            # match previous particles to current

            self.match_previous_particle_to_current(snapshot)

        # print("particleData: ",self.particleData)

        self.current_snapShotIndex += 1

    def correct_learning_rate(self, learning_rates_data: list):
        for i in range(len(learning_rates_data)):
            self.corrected_learning_rates_data[i] = learning_rates_data[i]
        self.learning_rate_corrected = True
 
    def find_array_in_list(self,target, list_of_arrays):
        for idx, arr in enumerate(list_of_arrays):
            # print("equal between: ",target, " and ", arr, " is: ",np.array_equal(target, arr))
            if np.array_equal(target, arr):

                return True
        return False

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

    def find_closest_particle(self, particle, shot, closest_rank=1):

        distances = np.linalg.norm(shot - particle, axis=1)

        minDistanceIndex = np.argmin(distances)

        if closest_rank > 1:
            for i in range(closest_rank - 1):
                minDistanceIndex = np.argmin(np.delete(distances, minDistanceIndex))
        return minDistanceIndex

    def get_default_learning_rate(self):
        return self.learning_rate_2D

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

    def get_original_shotData(self):
        return self.original_shotData

    def get_particle_id(self, particle,snapshotID, closest_rank=1,tolerance=0.01):
        target_snapshot = self.shotData[snapshotID]
        closest_particle_id_in_shot = self.find_closest_particle(particle, np.array(self.shotData[snapshotID]),closest_rank)

        # print("particle in get particle id: ",particle)
        # print("closest_particle_coor: ",closest_particle_id_in_shot)
        # print("target_snapshot: ",target_snapshot)
        # print("particleData with id: ",self.particleData[closest_particle_id_in_shot])
        # print("snapshotID relative: ",particle_relative_shotID)
        # print("particleData on this shot: ",self.particleData[closest_particle_id_in_shot]['coords'][particle_relative_shotID])
        for particle_id in self.particleData_2D:

            # print("iterating at particle_id: ",particle_id)
            # print("particleData_individual: ", self.particleData[particle_id]['coords'][snapshotID])
            # print("target_snapshot[closest_particle_id_in_shot]: ",target_snapshot[closest_particle_id_in_shot])
            # print("self.particleData[particle_id]['coords'][snapshotID]: ",self.particleData[particle_id]['coords'])

            particle_relative_shotID = self.find_relative_snapshotIndex(particle_id, snapshotID)
            distance = np.linalg.norm(target_snapshot[closest_particle_id_in_shot] - self.particleData_2D[particle_id]['coords'][particle_relative_shotID])
            # print("particle_relative_shotID: ",particle_relative_shotID)    
            if distance < tolerance:
                # print("found the particle id: ",particle_id)
                # particle_relative_shotID = self.find_relative_snapshotIndex(particle_id, snapshotID)
                # print("particle_relative_shotID: ",particle_relative_shotID)
                print(self.particleData_2D[particle_id]['coords'][particle_relative_shotID])
            
                return particle_id, self.particleData_2D[particle_id]['coords'][particle_relative_shotID]
            

                
        print("not found")
        print("orignal particle: ", target_snapshot[closest_particle_id_in_shot], " minus: ", self.particleData_2D['coords'][-1])
        KeyError("particle_id not found")

    def get_particle_id_from_available_ids(self,particle, snapshotID,id_list):
        target_coordinates = []
        
        for particle_id in id_list:
            particle_relative_shotID = self.find_relative_snapshotIndex(particle_id, snapshotID)
            target_coordinates.append(self.particleData_2D[particle_id]['coords'][particle_relative_shotID])
        closest_particle_id_in_shot = id_list[self.find_closest_particle(particle, np.array(target_coordinates))]
        return closest_particle_id_in_shot, self.particleData_2D[closest_particle_id_in_shot]['coords'][particle_relative_shotID]

    def get_particle_id_from_unmatched_ids(self,particle, snapshotID,matched_id_list):
        id_list = list(range(0, len(self.particleData_2D)))
        for id in matched_id_list:
            id_list.remove(id)
        target_coordinates = []
        for particle_id in id_list:
            particle_relative_shotID = self.find_relative_snapshotIndex(particle_id, snapshotID)
            # print("particle_relative_shotID: ",particle_relative_shotID)
            target_coordinates.append(self.particleData_2D[particle_id]['coords'][particle_relative_shotID])
        # print("target_coordinates: ",target_coordinates)
        closest_particle_id_in_shot = id_list[self.find_closest_particle(particle, np.array(target_coordinates))]
        
        return closest_particle_id_in_shot, self.particleData_2D[closest_particle_id_in_shot]['coords'][particle_relative_shotID]

    def get_total_num_of_particles(self):
        return len(self.particleData_2D)

    def get_particle_data(self):
        return self.particleData_2D

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

    def match_previous_particle_to_current(self, current_shot):

        def is_motion_random(historical_vel, observed_vel, motion_randomness):
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
            print("closest_particles: ",closest_particles)
            
            previous_index = closest_particles[3]
            current_index = closest_particles[2]
            # previous_index = closest_particles[2]
            # current_index = closest_particles[1]
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
          
            # explicit loop to remove the element from the current shot list
            for idx, particle in enumerate(current_shot_remain):
                if np.array_equal(particle, current_particle_to_match):
                    del current_shot_remain[idx]
                    break
            
            print("we are matching: ",prev_particle_coor, " with ", current_particle_to_match, " particle id: ",previous_particle_id)


            if previous_particle_id not in self.particleData_2D:
                self.particleData_2D[previous_particle_id] = {'coords': [], 'snapshotIndexList': [], 'snapshotIndexSet': set(), 'learning_rate': []}

            


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

                        # print()

                        previous_2D_shots_selected = self.particleData_2D[previous_particle_id]['coords'][previous_relative_index - self.NOS_per_section : previous_relative_index]
                        print("length of previous_2D_shots_selected: ",len(previous_2D_shots_selected))

                        # 
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
                        print("not enough data to do the reconstruction or feature disabled, use 2D historical velocity instead")
                        # calculate position from learning factor
                        if self.learning_rate_corrected:
                            learning_rates_2D = self.corrected_learning_rates_data[previous_particle_id]
                            print("learning_rates_2D: ",learning_rates_2D)
                            learning_rate_2D = learning_rates_2D[self.current_snapShotIndex]
                        else:
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

    def particle_projection(self, alpha, r_0):
        r_0_rotated=rotation(r_0,alpha)
        _, _, _, _, _, SOD, ODD,_,_ = self.reconstruction_conditions
        M_p = (SOD + ODD) / (SOD + r_0_rotated[1])
        
    
        return np.array([M_p * r_0_rotated[0], M_p * r_0_rotated[2]])

    def save_initial_particles(self, snapshot):
        for particle in snapshot:
            particle_id = self.assign_particle_id()
            # if self.paricleData is None:
            #     self.paricleData = {particle_id: particle}
            if particle_id not in self.particleData_2D:
                self.particleData_2D[particle_id] = {'coords': [], 'snapshotIndexList': [], 'snapshotIndexSet': set(), 'learning_rate':[]}
                
            self.particleData_2D[particle_id]['coords'].append(particle)
            self.particleData_2D[particle_id]['snapshotIndexList'].append(self.current_snapShotIndex)
            self.particleData_2D[particle_id]['snapshotIndexSet'].add(self.current_snapShotIndex)
            self.particleData_2D[particle_id]['learning_rate'].append(learning_rate_2D)

    def rank_particle_distances(self, previous_snapshot, current_snapshot, search_radius):
        # new method, ranked by the line distances
        # established a heap queue
        ranked_particle_heapq = []

        for i in range(len(previous_snapshot)):
            for j in range(len(current_snapshot)):
                equation_end_point = previous_snapshot[i]
                current_point = current_snapshot[j]

                # get the line of possible existance
                line_in_new_basis = self.particle_to_line_in_new_basis(equation_end_point)
                
                print("line_in_new_basis: ",line_in_new_basis)
                print("current_snapshot[j]: ",current_snapshot[j])
                # get the closest particle in the current snapshot to the line
                distance = self.distance_to_line(current_point, line_in_new_basis)

                # distance = np.linalg.norm(current_point - equation_end_point)

                # if the particles are more in the middle, we rely more on close neighbor than matrix transformation
                if abs(current_point[0]) < 3 or abs(current_point[1]) < 3:
                    learning_factor_stability = 0.9
                else:
                    learning_factor_stability = 0.6

                learning_factor_stability = 0.5/abs(current_point[0]) + 0.5/abs(current_point[1])

                # store the particles together by priority of its line distances. j is the index of the curret particle, i is the index of the previous particle
                heapq.heappush(ranked_particle_heapq, ((abs(current_point[1] - equation_end_point[1]) + abs(current_point[0] - equation_end_point[0]))* learning_factor_stability + abs(distance *(1-learning_factor_stability)) , abs(current_point[1] - equation_end_point[1]), j, i))
            


        # #  old method ranked by closest neighbor between particles
        # # turn previous snapshot into a N by 2 matrix
        # previous_snapshot = np.array(previous_snapshot).reshape(-1,2)

        # # grab the x axis
        # previous_x = previous_snapshot[:,0]
        # # grab the y axis
        # previous_y = previous_snapshot[:,1]

        # # turn current snapshot into a N by 2 matrix
        # current_snapshot = np.array(current_snapshot).reshape(-1,2)

        # # grab the x axis
        # current_x = current_snapshot[:,0]
        # # grab the y axis
        # current_y = current_snapshot[:,1]

        # # get the large matrix to fix the case where there are different number of particles in the previous and current snapshot
        # x_large_matrix = np.tile(previous_x, (len(current_x), 1)) - np.tile(current_x, (len(previous_x), 1)).T
        # y_large_matrix = np.tile(previous_y, (len(current_y), 1)) - np.tile(current_y, (len(previous_y), 1)).T

        # # print("x_large_matrix: ",x_large_matrix)
        # # print("y_large_matrix: ",y_large_matrix)

        # # get the distance matrix
        # distance_matrix = np.sqrt(x_large_matrix**2 + y_large_matrix**2)
        # # print("distance_matrix: ",distance_matrix)

        # # create a heap queue to store the ranked particles
        # ranked_particle_heapq = []

        # # Get the sorted indices of the flattened distance_matrix
        # sorted_indices = np.argsort(distance_matrix.ravel())

        # # Convert the flattened indices to 2D row and column indices
        # row_col_indices = np.unravel_index(sorted_indices, distance_matrix.shape)

        # # Print the values in distance_matrix in ascending order along with their row and column indices
        # for i in range(len(sorted_indices)):
        #     row, col = row_col_indices[0][i], row_col_indices[1][i]
        #     # print(f"Value: {distance_matrix[row, col]}, Row: {row}, Col: {col}")
        #     distance = distance_matrix[row, col]
        #     # store the closest particles together. row is the index of the curret particle, col is the index of the previous particle
        #     heapq.heappush(ranked_particle_heapq, (distance, 0,row, col))
        #     # print("original ranked_particle_heapq: ",len(ranked_particle_heapq))

        return ranked_particle_heapq

    def particle_to_line_in_new_basis(self,projected_particle_coodinates): 
        # magnification_factor = self.SOD / (self.SOD + self.ODD)
        source_3D_coorindates = [0, -self.SOD, 0] # sitting on the y axis in every shot
        projected_particle_3D_coodinates = [projected_particle_coodinates[0], self.RDD, projected_particle_coodinates[1]]

        # transform the particles to the new rotated basis
        transformed_source_3D_coorindates = rotation(source_3D_coorindates, -self.alpha) 
        transformed_projected_particle_3D_coodinates = rotation(projected_particle_3D_coodinates, -self.alpha)

        new_source_3D_coorindates = source_3D_coorindates
        # calculate the line equation end particles in the new basis

        projected_transformed_source_3D_coorindates = self.find_intersection(new_source_3D_coorindates,transformed_source_3D_coorindates, self.RDD)
        projected_transformed_projected_particle_3D_coodinates = self.find_intersection(new_source_3D_coorindates,transformed_projected_particle_3D_coodinates, self.RDD)
        # print("projected_transformed_source_3D_coorindates: ",projected_transformed_source_3D_coorindates)
        # print("projected_transformed_projected_particle_3D_coodinates: ",projected_transformed_projected_particle_3D_coodinates)
        transformed_source_2D_coorindates = np.array((projected_transformed_source_3D_coorindates[0],projected_transformed_source_3D_coorindates[2]))
        transformed_projected_particle_2D_coodinates = np.array((projected_transformed_projected_particle_3D_coodinates[0],projected_transformed_projected_particle_3D_coodinates[2]))



        return np.array([transformed_source_2D_coorindates, transformed_projected_particle_2D_coodinates])

    # Function to calculate distance from a point to the line
    def distance_to_line(self, point, ABMatrix):
        line_point1, line_point2 = ABMatrix
        line_point1 = np.array([float(elem) for elem in line_point1])
        line_point2 = np.array([float(elem) for elem in line_point2])
        line_vec = line_point2 - line_point1
        point_vec = point - line_point1
      
        # print("line_vec: ",line_vec)
        # print("Type of line_vec: ",type(line_vec))
        # print("type of lin vec element: ",type(line_vec[0]))
        
        line_len = np.linalg.norm(line_vec)
        line_unitvec = line_vec / line_len
        point_vec_scaled = point_vec / line_len
        t = np.dot(line_unitvec, point_vec_scaled)    
        nearest = line_point1 + t * line_vec
        return np.linalg.norm(nearest - point)
    
    def find_intersection(self, p1, p2, RDD):
        # Convert points to numpy arrays
        p1 = np.array(p1)
        p2 = np.array(p2)

        # Direction vector of the line
        direction = p2 - p1

        # Check if the line is parallel to the plane
        if direction[1] == 0:
            return "No intersection or line lies in the plane"

        # Calculate the t parameter for the line equation
        t = (RDD - p1[1]) / direction[1]

        # Calculate the intersection point
        intersection_point = p1 + t * direction

        return intersection_point

    




# In[ ]:


# 3D reconstruction functions

def generateEstimatedPositions(alpha, proj_used_index, N, xz_proj, conditions):
    _,delta_T, _, theta_degree, _, SOD, ODD,method,dataPiling = conditions
    theta = np.deg2rad(theta_degree)
    positions_predicted = np.zeros((N, 3))
    
    #     # Record the start time
    # start_time = time.time() 

    # use the new method, more time, more accuracy
    # values_this_round = proj2r0_acc(xz_proj[proj_used_index : proj_used_index+N-2, :], theta, SRD, RDD, delta_T)
    # position_rotated = position_rotated[0]
    # velocity_rotated = velocity_rotated[0]
    # acc_rotated = acc_rotated[0]

    # old method for time efficiency

    # print("proj_used_index: ",proj_used_index)
    # print(proj_used_index+N-2)
    # print("N: ",N)

    # print("xz_proj[proj_used_index : proj_used_index+N-2, :]: ",xz_proj[proj_used_index : proj_used_index+N-2, :])
    # print("proj_used_index: ",proj_used_index)
    # print("proj_used_index+N-2: ",proj_used_index+N-2)
    values_this_round = proj2r0_acc_old(xz_proj[proj_used_index : proj_used_index+N-2, :], theta, SOD, ODD, delta_T)
        # Record the end time
    # end_time = time.time()

    # # Calculate and print the total runtime
    # runtime = end_time - start_time
    # print(f"The runtime of projr20 is {runtime} seconds.")
    if method == 'acceleration':
        x0, y0, z0, u, v, w, a_x, a_y, a_z = values_this_round
        # print("values_this_round",values_this_round)
        position_rotated =  np.transpose(rotation([x0, y0, z0], alpha))
        # print("position_rotated",position_rotated)
        x0, y0, z0 = position_rotated[0][0], position_rotated[0][1], position_rotated[0][2]
        positions_predicted[0, :] = position_rotated
        velocity_rotated = np.transpose(rotation([u, v, w], alpha))
        u, v, w = velocity_rotated[0][0], velocity_rotated[0][1], velocity_rotated[0][2]
        acc_rotated =  np.transpose(rotation([a_x, a_y, a_z], alpha))
        a_x, a_y, a_z = acc_rotated[0][0], acc_rotated[0][1], acc_rotated[0][2]


        for j in range(1, N):
            time = delta_T * (j)
            positions_predicted[j, :] = [x0+u*time+0.5*a_x*time**2, y0+v*time+0.5*a_y*time**2, z0+w*time+0.5*a_z*time**2]


    elif method == 'linear':
        x0, y0, z0, u, v, w,_,_,_ = values_this_round
        position_rotated = np.transpose(rotation([x0, y0, z0], alpha))
        x0, y0, z0 = position_rotated
        positions_predicted[0, :] = position_rotated
        velocity_rotated = np.transpose(rotation([u, v, w], alpha))
        u, v, w = velocity_rotated

        for j in range(1, N):
            time = delta_T * (j)
            positions_predicted[j, :] = [x0+u*time, y0+v*time, z0+w*time]

    return positions_predicted   

def proj2r0_acc_old(proj, theta, SOD, ODD, delta_T):
    NOS = len(proj)
    print("proj: ",proj)
    SDD = SOD + ODD
    row_number_A = 2 * NOS + 2 * (NOS - 1)
    col_number_A = 1 + 2 * NOS
    # print("row_number_A: ",row_number_A)
    # print("col_number_A: ",col_number_A)
    A = np.zeros((row_number_A, col_number_A))
    b = np.zeros((row_number_A,1))
    
    for j in range(NOS):
        xi_j, zi_j = proj[j]
        A[2*j, 0] = 1
        A[2*j, 2*j+2] = -zi_j / SDD
        b[2*j] = zi_j * SOD / SDD
        
        A[2*j+1, 2*j+1] = -1
        A[2*j+1, 2*j+2] = xi_j / SDD
        b[2*j+1] = -xi_j * SOD / SDD
    
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
    r0 =[x[1], x[2], x[0], x[new_col_num-6], x[new_col_num-5], x[new_col_num-4], x[new_col_num-3], x[new_col_num-2], x[new_col_num-1]]
    return r0

def proj2r0_acc(xz_proj, theta, SOD, ODD, dt):
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

def Phase4_trace_3d(conditions, xz_proj):
    _,delta_T, NOS, theta_degree, N, SOD, ODD,method,dataPiling = conditions

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
            # print("project Used indes after: ",proj_used_index)
            alpha = -theta*(proj_used_index)
            # print("I ran it again after")
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
                # print("while project used index is less than NOS larger loop: ", xz_proj)
                new_positions = generateEstimatedPositions(alpha, proj_used_index, NOS_per_section,xz_proj,conditions)
                # print('prev',prev_proj_index)
                # print('proj_index',proj_used_index)
                # print('last positions',last_positions)
                # print('new positions',new_positions)
                combined_positions = np.concatenate([(new_positions[:len(last_positions), :] + last_positions)/2, new_positions[len(last_positions):, :]], axis=0)

                positions_predicted[proj_used_index : proj_used_index+NOS_per_section, :] = combined_positions
                # print("retro called")
                proj_used_index += NOS_per_section + 1
                # print("proj_used_index_new: ",proj_used_index)
                

    elif dataPiling == 'overlap':
        for i in range(round(NOS - N)):
            print("iterations: ",i)
            alpha = -theta * (proj_used_index - 1)  # alpha is for tracking the degree rotated from the 1st shot

            if proj_used_index == 1:
                positions_predicted[proj_used_index:N+1, :] = generateEstimatedPositions(alpha, proj_used_index, NOS_per_section,xz_proj,conditions)
            else:
                # take every N shots from every index, and take average of them
                last_positions = positions_predicted[proj_used_index: proj_used_index + prev_NOS_section - 1, :]

                
                new_positions = generateEstimatedPositions(alpha, proj_used_index, NOS_per_section,xz_proj,conditions)
        


                new_positions[:len(last_positions), :] = (new_positions[:len(last_positions), :] + last_positions) / 2
                positions_predicted[proj_used_index:proj_used_index + new_positions.shape[0], :] = new_positions
      

            proj_used_index += 1


    return positions_predicted

def smooth_points(estimated_positions, method, frame_size):
    """
    Smoothens the given 3D estimated positions using one of the specified methods.
    
    Parameters:
    - estimated_positions (ndarray): Nx3 array of estimated 3D positions.
    - method (str): Smoothing method ('avg', 'sg', or 'cb').
    - frame_size (int): Window size for moving average or Savitzky-Golay filter.
    
    Returns:
    - ndarray: Nx3 array of smoothened 3D positions.
    """
    # Initialize the filtered measurements with the original data
    filtered_measurements = estimated_positions.copy()
    
    # Moving Average
    if method == 'avg':
        for i in range(3):  # Loop over each dimension
            filtered_measurements[:, i] = np.convolve(estimated_positions[:, i], np.ones(frame_size)/frame_size, mode='same')
    
    # Savitzky-Golay Filter
    elif method == 'sg':
        for i in range(3):  # Loop over each dimension
            filtered_measurements[:, i] = savgol_filter(estimated_positions[:, i], frame_size, 2)
    
    # Cubic Smoothing Spline
    elif method == 'cb':
        x = np.arange(estimated_positions.shape[0])
        for i in range(3):  # Loop over each dimension
            spl = UnivariateSpline(x, estimated_positions[:, i], s=0.5)
            filtered_measurements[:, i] = spl(x)
    
    return filtered_measurements


# In[ ]:


# graphing functions
def number_to_binary_list(number):
    binary_str = bin(number)[2:]  # Convert to binary and remove the '0b' prefix
    binary_str = binary_str.zfill(3)  # Pad with zeros to make sure it has 3 digits
    binary_list = np.array([int(b) for b in binary_str] ) # Convert each binary digit to integer
    return binary_list

# def plotting_single(positions, particle_id,ax, label='particle'):
#     NOS, _ = positions.shape
#     particle_id=int(particle_id)+1
#     if particle_id<7 and particle_id>=0:
#         col=number_to_binary_list(int(particle_id))
#     else:
#         col=np.array([random.uniform(0.6, 1) for _ in range(3)])
#         print(col)

    
#     # Plot the curve with gradually changing color
#     for i in range(NOS - 1):
#         color_brightness = (i/(NOS-1))*col
        
#         x1, y1, z1 = positions[i]
#         x2, y2, z2 = positions[i+1]
        
#         # Draw a line segment with the computed color
#         ax.plot([x1, x2], [y1, y2], [z1, z2], color=color_brightness)
#     ax.plot([x1, x2], [y1, y2], [z1, z2], color=color_brightness,label=label+str(particle_id-1))
#     # Label the axes
#     ax.set_xlabel('X')
#     ax.set_ylabel('Y')
#     ax.set_zlabel('Z')
#     plt.legend()

def plotting_single(fig, positions, particle_id):
    # print("positions: ",positions)
    NOS, _ = positions.shape
    particle_id = int(particle_id) + 1
    if particle_id < 7 and particle_id >= 0:
        col = number_to_binary_list(int(particle_id))
    else:
        col = np.array([random.uniform(0.6, 1) for _ in range(3)])

    # # Create a Plotly figure
    # fig = go.Figure()

    # Plot the curve with gradually changing color
    # for i in range(NOS - 1):
        # color_brightness = (i / (NOS - 1)) * col
        # x1, y1, z1 = positions[i]
        # x2, y2, z2 = positions[i + 1]


        # Add a line segment to the figure
    fig.add_trace(go.Scatter3d(
        x=positions[:, 0],
        y=positions[:, 1],
        z=positions[:, 2],
        mode='lines',
        name='Particle ' + str(particle_id - 1),
        # showlegend=True  # Only for traces you want in the legend
        # line=dict(color='rgb({},{},{})'.format(*color_brightness*255),
        #           width=10)
    ))

    # Update layout for axes labels and legend
    fig.update_layout(
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z'
        ),
        legend_title_text='Particle ' + str(particle_id - 1)
    )
    
    return fig


# def show_in_window(fig):
#     import sys, os
#     import plotly.offline
#     from PyQt5.QtCore import QUrl
#     from PyQt5.QtWebEngineWidgets import QWebEngineView
#     from PyQt5.QtWidgets import QApplication
    
#     plotly.offline.plot(fig, filename='name.html', auto_open=False)
    
#     app = QApplication(sys.argv)
#     web = QWebEngineView()
#     file_path = os.path.abspath(os.path.join(os.path.dirname("Phase7_Python"), "name.html"))
#     web.load(QUrl.fromLocalFile(file_path))
#     web.show()
#     sys.exit(app.exec_())


# In[ ]:


# file processing function
def rename_files_replace_space(directory_path):
    # Get a list of all files in the directory
    filenames = os.listdir(directory_path)
    
    for filename in filenames:
        # Replace spaces with underscores
        new_filename = filename.replace('_', '')
        
        # Construct the full old and new file paths
        old_filepath = os.path.join(directory_path, filename)
        new_filepath = os.path.join(directory_path, new_filename)
        
        # Rename the file
        os.rename(old_filepath, new_filepath)


# 

# In[ ]:


# main function starts here
# define different learning rates for investigating different sorting and smoothing model

alpha = np.radians(theta_degrees)  # Example rotation angle in radians

rates_conditions = [learning_rate_2D, motion_randomness, learning_rate_3D, NOS_per_section]
path_finder = pf(alpha,rates_conditions,conditions)

print("list of files: ",os.listdir(folderName))
sorted_filenames = sorted(os.listdir(folderName), key=lambda x: int(x.split(file_prefix)[1].split('.csv')[0]))
# print(len(sorted_filenames[0]))
k = 0
for fileIndex in range(min(len(sorted_filenames), NOS)):
    file = sorted_filenames[fileIndex]
    print("file: ",file)
    if file.endswith(".csv"):
        filename = os.path.join(folderName, file)
        input_data = pd.read_csv(filename, header=None)
        input_data = np.array(np.transpose(input_data))
        values =  input_data[0]
        print("values:", values)
        print("read file: ", filename)

        paired_values = []
        i = 0
        for j in range(len(values)//2):
            
            inputList = (values[i:i+2]- offset)*pixelResolution
           
            # input format, list of tuple of two elements (x,y)
            paired_values.append(inputList)

            # scambled_values = random.shuffle(paired_values.copy())
            # print("paired values:", paired_values)
            
            i += 2

        # print("paired values:", paired_values)
        # paired_values.append(xz_proj[k])
        path_finder.append(paired_values)
        k += 1

# Assuming path_finder.get_particle_data() returns your data as a dictionary
sorted_particle_data = path_finder.get_particle_data()
print("sorted_particle_data: ",sorted_particle_data)

NumOfDataPoints = len(sorted_particle_data)
print("NumOfDataPoints: ", NumOfDataPoints)

print(sorted_particle_data)
estimated_positions = np.zeros((NOS,3*NumOfDataPoints))


for i in range(NumOfDataPoints):
    # print("inline, ",np.array(sorted_particle_data[i]['coords']))
    estimated_positions_single = Phase4_trace_3d(conditions, np.array(sorted_particle_data[i]['coords']))
    # print("completed")
    estimated_positions_single = smooth_points(estimated_positions_single, 'sg',NOS_per_section)

    estimated_positions[:,i*3:i*3+3] = estimated_positions_single

estimated_positions_byPass = np.zeros((NOS,3*NumOfDataPoints))
# extra synthetic data points (bypassing sorting)
for i in range(NumOfDataPoints):
    estimated_positions_single = smooth_points(Phase4_trace_3d(conditions, xz_proj[:,i*2:i*2+2]), 'sg',NOS_per_section)
    print("count",i)
    estimated_positions_byPass[:,i*3:i*3+3] = estimated_positions_single

estimated_positions_new = estimated_positions
estimated_positions_graph = np.column_stack((estimated_positions, estimated_positions_byPass, real_positions))



# In[ ]:


# %% recalculate learning rates
learning_rate_floor = path_finder.get_default_learning_rate()
learning_rate_floor = 1
learning_rate_ceiling = 1

# empirical, adjust here
min_distance_to_origin = 0.6
max_distance_to_origin = 6

learning_rate_lists = {}

for i in range(path_finder.get_total_num_of_particles()):
    estimated_positions_individual = estimated_positions_new[:,i*3:i*3+3]

    learning_rate_lists[i] = map_range(np.linalg.norm(estimated_positions_individual, axis=1), min_distance_to_origin, max_distance_to_origin, learning_rate_floor, learning_rate_ceiling)

print("learning_rate_lists: ",learning_rate_lists)

corrected_path_finder = pf(alpha,rates_conditions,conditions)

corrected_path_finder.correct_learning_rate(learning_rate_lists)

shotData = path_finder.get_original_shotData()
for i in range(len(shotData)):
        # print("shotData",self.shotData[i])
        corrected_path_finder.append(shotData[i])
# corrected_path_finder.re_run(path_finder.get_original_shotData())

sorted_particle_data_corrected = corrected_path_finder.get_particle_data()

# Convert the NumPy matrix to a Pandas DataFrame
outputMatrix = sorted_particle_data_corrected[0]['coords']
for i in range(1,NumOfDataPoints):
     outputMatrix = np.column_stack((outputMatrix, sorted_particle_data_corrected[i]['coords']))
df = pd.DataFrame(outputMatrix)
# Export the DataFrame to an Excel file
df.to_excel("sorted_particle_data_corrected.xlsx", index=False)
# Convert the NumPy matrix to a Pandas DataFrame
df = pd.DataFrame(xz_proj)
# Export the DataFrame to an Excel file
df.to_excel("xz_proj.xlsx", index=False)

estimated_positions_corrected = np.zeros((NOS,3*NumOfDataPoints))
for i in range(NumOfDataPoints):
    print("sorted_particle_data_corrected[i]['coords']",sorted_particle_data_corrected[i]['coords'])
    estimated_positions_single = Phase4_trace_3d(conditions, np.array(sorted_particle_data_corrected[i]['coords']))
    estimated_positions_single = smooth_points(estimated_positions_single, 'sg',NOS_per_section)

    estimated_positions_corrected[:,i*3:i*3+3] = estimated_positions_single

# estimated_positions_corrected_new = estimated_positions_corrected

# estimated_positions_corrected_graph = np.column_stack((estimated_positions_corrected,estimated_positions_byPass, real_positions))

estimated_positions_corrected_graph = np.column_stack((estimated_positions_corrected,real_positions))
# estimated_positions_corrected_graph= estimated_positions_corrected_new

print("estimated_positions_corrected_graph: ",estimated_positions_corrected_graph)

# Convert the NumPy matrix to a Pandas DataFrame
df = pd.DataFrame(estimated_positions_corrected_graph)

# Export the DataFrame to an Excel file
# df.to_excel("estimated_positions_corrected_graph.xlsx", index=False)

# NumOfDataPoints += 2


# In[ ]:


print("shape of estimated_positions_graph: ",estimated_positions_graph.shape)
import plotly.graph_objs as go
# Create the figure and axes
# fig = plt.figure()
fig = go.Figure()
# fig.update_layout(
#     scene=dict(
#         xaxis_title='X',
#         yaxis_title='Y',
#         zaxis_title='Z'
#     ),
#     legend_title_text='Particles'
# )
# ax = fig.add_subplot(111, projection='3d')
_, column = estimated_positions_graph.shape
for i in range(column//3):
    # plotting_single(estimated_positions_graph[:,i*3:i*3+3],i,ax)
    plotting_single(fig, estimated_positions_graph[:,i*3:i*3+3],i)

# NumOfDataPoints -=2
# print("NumOfDataPoints: ",NumOfDataPoints)

fig2 = go.Figure()
# ax2 = fig2.add_subplot(111, projection='3d')

x_min = y_min = z_min = float('inf')
x_max = y_max = z_max = float('-inf')
print(NumOfDataPoints)
_, column = estimated_positions_corrected_graph.shape
print("estiamted corrected shape:",estimated_positions_corrected_graph.shape)
for i in range(column//3):
    plotting_single(fig2, estimated_positions_corrected_graph[:,i*3:i*3+3],i)
    x_min = min(min(estimated_positions_corrected_graph[:,i*3]), x_min)
    x_max = max(max(estimated_positions_corrected_graph[:,i*3]), x_max)
    y_min = min(min(estimated_positions_corrected_graph[:,i*3+1]), y_min)
    y_max = max(max(estimated_positions_corrected_graph[:,i*3+1]), y_max)
    z_min = min(min(estimated_positions_corrected_graph[:,i*3+2]), z_min)
    z_max = max(max(estimated_positions_corrected_graph[:,i*3+2]), z_max)

    # print("estimated_positions_corrected_graph[:,i*3:i*3+3]",estimated_positions_corrected_graph[:,i*3:i*3+3])
# plotting_single(estimated_positions_corrected_graph,i,ax2)

print(x_min, x_max, y_min, y_max, z_min, z_max)
# if x_max - x_min < 1e-3: 
#     ax2.set_xlim3d(x_min-1, x_max+1)

# if y_max - y_min < 1e-3:
#     ax2.set_ylim3d(y_min-1, y_max+1)

# if z_max - z_min < 1e-3:
#     ax2.set_zlim3d(z_min-1, z_max+1)

# show_in_window(fig)
# show_in_window(fig2)
fig.show()
fig2.show()


