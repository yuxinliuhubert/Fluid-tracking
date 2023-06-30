%% Data Input
function [real_positions,positions_predicted] = Phase4_pt_3d(initial_position_3d, noise, delta_T, NOS,theta_degree,N)
%N means the number of shots to use for one matrix; NOS/N must be a whole
%number
format default
SRD = 1; % m, Source-Reference Distance
RDD = 1; % m, Reference-Detector (screen) Distance
% % theta_degree = 10; % clock-wise degree, what is the angle of camera rotation before each shot
% % NOS = 100; % number of shots
% % delta_T = time between shots
% % initial_position_3d = [x0;y0;z0], a column vector
% noise=5e-3; %the 68% chance deviation of the projection measurement from the perfect measurement

method = 0; % 0 for least square, 1 for kalman

v = @(t) [3*t+7;2;1*t+2.2]; %velocity function

% Process error
delta_P_X = 1e-3; % m
delta_P_VX = 1e-3; % m/sec

% Observation error
delta_O_X= 2e-3; % m
delta_O_VX = 1.5e-4; % m
% delta_O_Y = delta_O_X;
% delta_O_VY = delta_O_VX;

%% Calculation part 1: projection generation
theta = theta_degree/180*pi;
% % First shot: without rotation
r0_0=initial_position_3d;
real_positions = initial_position_3d;
M_p = (SRD+RDD)/(SRD+r0_0(2));
x_proj=M_p*r0_0(1)+randn*noise/2;
z_proj=M_p*r0_0(3)+randn*noise/2;
% % Second shot: with theta degrees rotation each time
for k = 1:NOS-1
    r0_k=r0_0+integral(v,0,delta_T*k,'ArrayValued', true); % true location in the original frame of reference
    real_positions = [real_positions; r0_0 + integral(v,0,delta_T*k,'ArrayValued', true)'];
    r_now=T(r0_k,theta*k);
    M_p = (SRD+RDD)/(SRD+r_now(2)); % magnification of particle
    x_proj=[x_proj;M_p*r_now(1)+randn*noise/2];
    z_proj=[z_proj;M_p*r_now(3)+randn*noise/2];
end
% x_proj and z_proj are NOS by 1 vectors with all the x-coordinates projection on the screen
xz_proj=[x_proj, z_proj];

proj_used_index=1;
positions_predicted=[];
%% Calculation part2: geting the measured values
for i=1:round(NOS/N)
    alpha=-theta*(proj_used_index-1);%alpha is for tracking the degree rotated from the 1st shot
    values_this_round=proj2r0_acc(xz_proj(proj_used_index:(proj_used_index+N-1),:),theta,SRD,RDD,delta_T);
    [x0, y0, z0 ,u ,v ,w ,a_x, a_y ,a_z]=deal(values_this_round(1),values_this_round(2),values_this_round(3),values_this_round(4),values_this_round(5),values_this_round(6),values_this_round(7),values_this_round(8),values_this_round(9));
    position_rotated=T([x0;y0;z0],alpha)';
    [x0, y0, z0]=deal(position_rotated(1), position_rotated(2) , position_rotated(3) );
    positions_predicted=[positions_predicted;position_rotated];
    velocity_rotated=T([u;v;w],alpha)';
    [u ,v ,w]=deal(velocity_rotated(1), velocity_rotated(2) , velocity_rotated(3) );
    acc_rotated=T([a_x; a_y ;a_z],alpha)';
    [a_x, a_y ,a_z]=deal(acc_rotated(1), acc_rotated(2) , acc_rotated(3) );
    for j =1:N-1
        time=delta_T*j;
        positions_predicted=[positions_predicted;x0+u*time+0.5*a_x*time^2, y0+v*time+0.5*a_y*time^2, z0+w*time+0.5*a_z*time^2 ];
    end

    proj_used_index=proj_used_index+N;
end

%% function for rotating coordinate axes by a given angle alpha
    function [r2]=T(r1,alpha) %r1=[x1;y1;z1]and r2=[x2;y2;z2] describe the same physical pt but in 2 coordinate systems with
        % theta degress rotation CCW
        r2=[cos(-alpha) -sin(-alpha) 0; sin(-alpha) cos(-alpha) 0;0 0 1]*r1;
        % note that we are rotating the coordinate axes, not the location vector of the particle
    end


%row_number_A=row_number_A+NOS ;%which is different from phase 1. We get a new equation each shot for the change of z coordinate due to velocity





end