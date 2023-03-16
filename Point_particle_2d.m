%% Point Particle Static Position Tracking in 2D (view from top)
% We assume that the particle is a point particle in 2D with almost zero velocity.  
%% Data Input
SRD = 1; % m, Source-Reference Distance
RDD = 1; % m, Reference-Detector (screen) Distance
theta_degree = 30; % clock-wise degree, what is the angle of camera rotation before each shot
NOS = 2; % number of shots

% True location of particle; could be negative number
x0 = 0.3; % m, relative to reference
y0 = 0.2; % m
% M_particle = (SRD+RDD)/(SRD+y0); % magnification of particle

% projected x input; could be used to create kalman filter with rotation
% matrix to minimize the effect of noise. This version assumes no noise
% from measurement
x_proj = [0.4; 0.5];

%% Calculations
% Generate vectors and constants
M = (SRD+RDD)/SRD; % magnification of reference
theta = theta_degree/180*pi;

% Matrix of transformation
A = [cos(theta) -sin(theta); sin(theta) cos(theta)];

for i = 1:NOS % process each shot
    % take a guess
    




end



%% Final Processing & Saving





%% Graphing

