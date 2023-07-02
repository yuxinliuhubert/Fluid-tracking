initial_positions = [0,0,0];
noise = 1e-5;
delta_T = 0.007;
NOS = 240;
theta_degrees = 1.8;
NOS_per_section = 30; % must be larger than 5 to satisfy equations
% v = @(t) [3*t+7;2;1*t+2.2]; %velocity function

v = @(t) (t<delta_T*round(NOS/2)).* [2;1*t+3*t^2;2] + ...
(t>=delta_T*round(NOS/2) & t < delta_T*round(NOS*0.8)) .* [3*t+7;2;1*t+2.2] + ...
(t>=delta_T*round(NOS*0.8)).* [2-4*t^2;1*t;2];

conditions = [noise, delta_T, NOS,theta_degrees,NOS_per_section];

[real_positions, estimated_positions] = Phase4_pt_3d(initial_positions,conditions,v)


Phase4Graph(real_positions, estimated_positions,conditions,v)

Phase3_pt_3d_success(initial_positions,noise,delta_T,NOS,theta_degrees, v)