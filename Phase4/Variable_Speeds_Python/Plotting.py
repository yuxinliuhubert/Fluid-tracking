import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plotting_single(positions, particle_id,ax):
    NOS, _ = positions.shape
    RGB=int(particle_id%3)
    # Plot the curve with gradually changing color
    for i in range(NOS - 1):
        col = np.full(3,i/(NOS-1))
        dummy=np.array([0,0,0])
        dummy[RGB]=1
        col=col*dummy
        x1, y1, z1 = positions[i]
        x2, y2, z2 = positions[i+1]
        
        # Draw a line segment with the computed color
        ax.plot([x1, x2], [y1, y2], [z1, z2], color=col)
    ax.plot([x1, x2], [y1, y2], [z1, z2], color=col,label='particle'+str(particle_id))
    # Label the axes
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.legend()
    

# # Create the figure and axes

# # Test the function with different data

# positions1 = np.array([[0.3993, 3.8213, 0.0], [2.0606, 3.9602, 0.1471], [3.3272, 3.3407, 0.2941], [3.9567, 2.0814, 0.4412], [3.8285, 0.4236, 0.5882], [2.9671, -1.3154, 0.7353], [1.5375, -2.8025, 0.8824], [-0.1865, -3.7529, 1.0294], [-1.8748, -3.9846, 1.1765], [-3.2041, -3.4533, 1.3235], [-3.9198, -2.2606, 1.4706], [-3.8848, -0.6351, 1.6176], [-3.1059, 1.1121, 1.7647], [-1.7323, 2.6463, 1.9118], [-0.0269, 3.6737, 2.0588], [1.6837, 3.9976, 2.2059], [3.0718, 3.5560, 2.3529], [3.8717, 2.4334, 2.5000], [3.9301, 0.8448, 2.6471], [3.2359, -0.9056, 2.7941], [1.9221, -2.4825, 2.9412], [0.2401, -3.5841, 3.0882], [-1.4878, -3.9993, 3.2353], [-2.9308, -3.6486, 3.3824], [-3.8126, -2.5993, 3.5294], [-3.9642, -1.0521, 3.6765], [-3.3567, 0.6965, 3.8235], [-2.1064, 2.3117, 3.9706], [-0.4527, 3.4843, 4.1176], [1.2876, 3.9896, 4.2647]])
# positions2 = np.array([[0.1997, 1.9107, 0.0], [1.0303, 1.9801, 0.1471], [1.6636, 1.6703, 0.2941], [1.9784, 1.0407, 0.4412], [1.9142, 0.2118, 0.5882], [1.4835, -0.6577, 0.7353], [0.7688, -1.4012, 0.8824], [-0.0933, -1.8764, 1.0294], [-0.9374, -1.9923, 1.1765], [-1.6020, -1.7266, 1.3235], [-1.9599, -1.1303, 1.4706], [-1.9424, -0.3175, 1.6176], [-1.5530, 0.5560, 1.7647], [-0.8661, 1.3231, 1.9118], [-0.0134, 1.8369, 2.0588], [0.8418, 1.9988, 2.2059], [1.5359, 1.7780, 2.3529], [1.9358, 1.2167, 2.5000], [1.9651, 0.4224, 2.6471], [1.6180, -0.4528, 2.7941], [0.9610, -1.2413, 2.9412], [0.1201, -1.7920, 3.0882], [-0.7439, -1.9996, 3.2353], [-1.4654, -1.8243, 3.3824], [-1.9063, -1.2996, 3.5294], [-1.9821, -0.5261, 3.6765], [-1.6784, 0.3482, 3.8235], [-1.0532, 1.1559, 3.9706], [-0.2264, 1.7421, 4.1176], [0.6438, 1.9948, 4.2647]])

# plotting_single(positions1,0, ax)
# plotting_single(positions2, 1,ax)


# # Show the plot
# plt.show()
