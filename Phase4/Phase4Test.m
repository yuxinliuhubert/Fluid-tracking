initial_positions = [0,0,0];
noise = 1e-5;
delta_T = 0.02;
NOS = 150;
theta_degrees = 1.8;
NOS_per_section = 7; % must be larger than 5 to satisfy equations
v = @(t) [3*t+7;2;1*t+2.2]; %velocity function
conditions = [noise, delta_T, NOS,theta_degrees,NOS_per_section];

[real_positions, estimated_positions] = Phase4_pt_3d(initial_positions,conditions,v)

Phase4Graph(real_positions, estimated_positions,conditions,v)