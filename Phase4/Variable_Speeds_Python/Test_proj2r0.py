import numpy as np
from proj2r0_acc_combination import proj2r0_acc_combination
x_proj=np.array([[0.9524], [3.2454], [-76.8505], [-2.9333], [-0.8678], [0.2900] ])
z_proj=np.array([[2.8571], [4.6202], [-87.4060], [-3.8000], [-2.2091], [-1.9513]])
xz_proj = np.column_stack((x_proj, z_proj))
print(np.math.pi)
print(proj2r0_acc_combination(xz_proj,np.math.pi/6,1,1,1))
