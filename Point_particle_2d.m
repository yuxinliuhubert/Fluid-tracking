%% Point Particle Static Position Tracking in 2D (view from top)
% We assume that the particle is a point particle in 2D with almost zero velocity.  
format long
%% Data Input
SRD = 1; % m, Source-Reference Distance
RDD = 1; % m, Reference-Detector (screen) Distance
theta_degree = 0; % clock-wise degree, what is the angle of camera rotation before each shot
D = 1e-4; % m, diameter of the fluid element, assuming 3D sphere/2D circle
delta_T = 1; % s, time between shots
a = [0 0]; % m/s^2, fluid element acceleration


% Process error
delta_P_X = 1e-3; % m
delta_P_VX = 1e-3; % m/sec

% Observation error
delta_O_X= 2e-3; % m
delta_O_VX = 1.5e-4; % m
delta_O_Y = delta_O_X;
delta_O_VY = delta_O_VX;

O_X = [delta_O_X^2, 0; 0, delta_O_VX^2];
O_Y = [delta_O_Y^2, 0; 0, delta_O_VY^2];

% % True location of particle; could be negative number
% x0 = 0.3; % m, relative to reference
% y0 = 0.2; % m
% M_particle = (SRD+RDD)/(SRD+y0); % magnification of particle

% projected x input; could be used to create kalman filter with rotation
% matrix to minimize the effect of noise. This version assumes no noise
% from measurement
% M >1
x_proj = [ % [particle left projection, particle right side projection]
    0.502-1.2*D/2, 0.502+1.2*D/2 ;
    0.502-1.2*D/2, 0.502+1.2*D/2 ;
    0.502-1.2*D/2, 0.502+1.2*D/2 ;
    0.502-1.2*D/2, 0.502+1.2*D/2 ;
    0.502-1.2*D/2, 0.502+1.2*D/2 ;
    0.502-1.2*D/2, 0.502+1.2*D/2 ;
    0.502-1.2*D/2, 0.502+1.2*D/2 ;
     0.502-1.2*D/2, 0.502+1.2*D/2 ;
    0.502-1.2*D/2, 0.502+1.2*D/2 ;
    0.502-1.2*D/2, 0.502+1.2*D/2 ;
    0.502-1.2*D/2, 0.502+1.2*D/2 ;
    0.502-1.2*D/2, 0.502+1.2*D/2 ;
    0.502-1.2*D/2, 0.502+1.2*D/2 ;
    0.502-1.2*D/2, 0.502+1.2*D/2 ;
    
    ]


%% Calculations
NOS = length(x_proj); % number of shots
% M = (SRD+RDD)/SRD; % magnification of reference

% Matrix of transformation; used to convert measurements from rotated basis
% to original basis. For every shot taken, A multiplies by once. 
theta = theta_degree/180*pi;
angle_transformation = [cos(theta) -sin(theta); sin(theta) cos(theta)];

% initialize x and y state vectors
y_state= [];
x_state = [];

% First shot
xmL = x_proj(1,1);
xmR = x_proj(1,2);
syms x0cL x0cR y0c
eqn = [(x0cL+D)/(y0c+SRD) == xmR/(RDD+SRD); x0cL/(y0c+SRD) == xmL/(SRD+RDD)];
S = solve(eqn,[x0cL, y0c]);
xpL = double(S.x0cL)
yp = double(S.y0c);
xp(1) = xpL + D/2;
v = [0 0]; % assume the fluid element isn't moving after first shot

delta_P_Y = delta_P_X * (SRD + RDD)/xmR; % calculating Y observation error, which derives from X observation error
% Take d/dt on both sides
delta_P_VY = delta_P_VX * (SRD + RDD)/xmR;

% Initial process covariance matrix 
P_X0 = [delta_P_X^2 delta_P_X*delta_P_VX; delta_P_X*delta_P_VX delta_P_VX^2];
P_Y0 = [delta_P_Y^2 delta_P_Y*delta_P_VY; delta_P_X*delta_P_VY delta_P_VY^2];

x_state(1,:) = [xp(1) v(1)]; % first index position, second process covariance
y_state(1,:) = [yp v(2)];

for i = 2:NOS % process each shot
    A = [1 delta_T; 0 1];
    B = [1/2*delta_T^2; delta_T]; % multiplied by acceleration 

    % take a guess
    x_kp = A * x_state(i-1,:)' + B*a(1)
    y_kp = A * y_state(i-1,:)' + B*a(2);

    % Calculate new covariance matrix
    P_X = A*P_X0*A' + 0;
    P_Y = A*P_Y0*A' + 0;

    % Calculate the kalman gain
    K_X = P_X*eye(2)/(eye(2)*P_X*eye(2) + O_X);
    K_Y = P_Y*eye(2)/(eye(2)*P_Y*eye(2) + O_Y);

    % Import Observation
    [e, vk] = getData(x_proj, i, D, SRD, RDD, angle_transformation, x_state, y_state, delta_T);
    xk = eye(2)*[e(1); vk(1)]
    xp(i) = e(1);
    yk = eye(2)*[e(2); vk(2)];

    % Caluclate current state
    x_state(i, :) = x_kp' + (K_X*(xk-eye(2)*x_kp))';
    y_state(i, :) = y_kp' + (K_Y*(yk-eye(2)*y_kp))';

    % Update process covariance martrix
    P_X0 = (eye(2) - K_X*eye(2))*P_X;
    P_Y0 = (eye(2) - K_Y*eye(2))*P_Y;



end

%% Final Processing & Saving
shot = (1:NOS)';
table = array2table([shot,xp', x_state, y_state],"VariableNames",["shot","x_proj","x_pos","v_x","y_pos","v_y"])



function [e,v] = getData(x_proj, index, D, SRD, RDD, angle_transformation, x_state, y_state, delta_T)
    xmL = x_proj(index,1);
    xmR = x_proj(index,2);
    syms x0cL y0c
    eqn = [(x0cL+D)/(y0c+SRD) == xmR/(RDD+SRD); x0cL/(y0c+SRD) == xmL/(SRD+RDD)];
    S = solve(eqn,[x0cL, y0c]);
    xpL = double(S.x0cL);

    e = angle_transformation^(index-1) * [xpL + D/2;double(S.y0c)];
    vx = (e(1) - x_state(index-1, 1))/delta_T;
    vy = (e(2) - y_state(index-1, 1))/delta_T;
    v = [vx vy];
 
end

%% Graphing

