close all
%% Input conditions
initial_positions = [0,0,0];
noise = 1e-4;
delta_T = 0.02;
NOS = 240;
theta_degrees = 38;
NOS_per_section = 8; % must be larger than 5 to satisfy equations
% v = @(t) [3*t+7;2;1*t+2.2]; %velocity function

v = @(t) (t<delta_T*round(NOS/2)).* [2;1*t+3*t^3;2] + ...
(t>=delta_T*round(NOS/2) & t < delta_T*round(NOS*0.8)) .* [3*t^3+7;2;1*t^2+2.2] + ...
(t>=delta_T*round(NOS*0.8)).* [2-4*t^2;1*t;2];

conditions = [noise, delta_T, NOS,theta_degrees,NOS_per_section];

%% Graphing Prep
% Get the size of the screen
screenSize = get(0, 'ScreenSize');
% Define the width and height for your figures
figWidth = screenSize(3) / 2;  % width is 1/3 of the screen width
figHeight = screenSize(4) / 2; % height is 1/2 of the screen height
% Create and position three figures side by side
f2 = figure('Position', [1, (screenSize(4)-figHeight)/2, figWidth, figHeight]);
f3 = figure('Position', [figWidth+1, (screenSize(4)-figHeight)/2, figWidth, figHeight]);


%% Automated Section
method = 'linear';
dataPiling = 'serial';
[real_positions, estimated_positions] = Phase4_pt_3d(initial_positions,conditions,v,method,dataPiling);
Phase4Graph(real_positions, estimated_positions,conditions,v,method,dataPiling,f2);


method = 'linear';
dataPiling = 'overlap';
[real_positions, estimated_positions] = Phase4_pt_3d(initial_positions,conditions,v,method,dataPiling);
Phase4Graph(real_positions, estimated_positions,conditions,v,method,dataPiling,f3);


% Phase3_pt_3d_success(initial_positions,noise,delta_T,NOS,theta_degrees, v);