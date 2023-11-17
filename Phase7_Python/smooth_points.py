
import numpy as np
from scipy.signal import savgol_filter
from scipy.interpolate import UnivariateSpline

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
