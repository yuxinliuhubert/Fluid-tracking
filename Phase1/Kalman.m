%% Data Input
function [kalman_position_velocity, r0_avg] = Kalman(true_position_3d, noise, delta_T, NOS,theta_degree, SRD, RDD)
format long
% SRD = 1; % m, Source-Reference Distance
% RDD = 1; % m, Reference-Detector (screen) Distance
% % theta_degree = 10; % clock-wise degree, what is the angle of camera rotation before each shot
% % NOS = 100; % number of shots
% % delta_T = 1; % s, time between shots
% % True location of particle; could be negative number
% x0 = 0.1; y0 = 0.1; z0=0.1;  % in m, relative to reference
x0 = true_position_3d(1);
y0 = true_position_3d(2);
z0 = true_position_3d(3);
% noise=5e-3; %the 68% chance deviation of the projection measurement from the perfect measurement 

a = [0 0 0]; % m/s^2, fluid element acceleration

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
r_now=[x0; y0;z0];
M_p = (SRD+RDD)/(SRD+r_now(2));
x_proj=[M_p*r_now(1)+randn*noise/2];
z_proj=[M_p*r_now(3)+randn*noise/2];
% % Second shot: with theta degrees rotation each time
for i = 1:NOS-1
    r_now=T(r_now,theta); % true location in the newly rotated axes
    M_p = (SRD+RDD)/(SRD+r_now(2)); % magnification of particle
    x_proj=[x_proj;M_p*r_now(1)+randn*noise/2];
    z_proj=[z_proj;M_p*r_now(3)+randn*noise/2];
end
% x_proj and z_proj are NOS by 1 vectors with all the x-coordinates projection on the screen
xz_proj=[x_proj, z_proj];


%% Calculation part2: getting the measured values 
xyz_measured=proj2r0(xz_proj,theta,SRD,RDD);

% initialize the states
xyz_initial = xyz_measured(1,:);
xyz_predicted(1,:) = xyz_initial;
velocity_predicted(1,:) = [0, 0, 0]; % x, y, z

% a state involves position and velocity. x_ki or x_kp will be a 2-by-1
% vector
% initial states
x_ki = [xyz_measured(1,1); 0];
y_ki = [xyz_measured(1,2); 0];
z_ki = [xyz_measured(1,3); 0];
    
% matrices set up 
A = [1, delta_T; 0, 1];
B = [0.5*delta_T^2; delta_T];

% initial process covariance matrix
P_ki = [delta_P_X^2, 0; 0, delta_P_VX^2];
% 2nd and 3rd entries are both assumed to be 0 instead of delta_P_X *
% delta_P_VX because the process errors of velocity and position are not
% really related to each other.

% observation error matrix R
R = [delta_O_X^2, 0; 0, delta_O_VX^2];

x_k_prev = x_ki;
y_k_prev = y_ki;
z_k_prev = z_ki;
P_k_prev = P_ki;

if (NOS > 2)
% iteration begins
for i = 2:height(xyz_measured) 
    
    % predicted states
    x_kp = A * x_k_prev + B * a(1);
    y_kp = A * y_k_prev + B * a(2);
    z_kp = A * z_k_prev + B * a(3);

    % predicted process cobariance matrix
    P_kp = A*P_k_prev*A' + 0; % 0 in place of Q
    P_kp(2:3) = 0;

    % calculating the kalman gain
    K = P_kp*eye(length(P_kp))/(eye(length(P_kp))*P_kp*eye(length(P_kp)) + R);

    % new observation
    x_O = [xyz_measured(i, 1); (xyz_measured(i,1) - xyz_measured(i-1,1))/delta_T] + 0; % + 0 in place of noise in the electronics
    y_O = [xyz_measured(i, 2); (xyz_measured(i,2) - xyz_measured(i-1,2))/delta_T];
    z_O = [xyz_measured(i, 3); (xyz_measured(i,3) - xyz_measured(i-1,3))/delta_T];

    % calculate current state
    x_k = x_kp + K*(x_O- eye(length(P_kp))*x_kp);
    y_k = y_kp + K*(y_O- eye(length(P_kp))*y_kp);
    z_k = z_kp + K*(z_O- eye(length(P_kp))*z_kp);

    % update the process covariance matrix
    P_k = (eye(length(P_kp)) - K*eye(length(P_kp)))*P_kp;

    % prepare the second round
    x_k_prev = x_k;
    y_k_prev = y_k;
    z_k_prev = z_k;
    P_k_prev = P_k;

end
else 
x_k = x_ki;
y_k = y_ki;
z_k = z_ki;

end


% Using the first measurement as a reference point, now run kalman filter
% on the subsequent measurements, assuming that the velocity and
% acceleration is 0 to start. For simplicity acceleration is set to 0
% the entire time. In other words, linear velocity is expected.


%% Part 3: System output
kalman_position_velocity = [x_k y_k z_k];
r0_avg=mean(xyz_measured);




 %% function for rotating coordinate axes by a given angle alpha 
   function [r2]=T(r1,alpha) %r1=[x1;y1]and r2=[x2;y2] describe the same physical pt but in 2 coordinate systems with
   % theta degress rotation CCW
   r2=[cos(-alpha) -sin(-alpha) 0; sin(-alpha) cos(-alpha) 0;0 0 1]*r1;
   % note that we are rotating the coordinate axes, not the location vector of the particle
   end
  
   %% 
   function [r0]= proj2r0(proj,theta,SRD,RDD) %takes in any N (even #) by 2 of projection coordinates and the angle 
   %theta btw the shots and outputs the predicted [x0, y0, z0] sets in coordinates when theta=0 (original
   %position before camera machine starts rotating)
   %r0=[x0_1, y0_1,z0_1; x0_2,y0_2,z_02;.....]

r0=zeros(height(proj)-1,3); %initialize the predicted positions
       for i = 1:(height(proj)-1)
           xi_1=proj(i,1); xi_2=proj(i+1,1); alpha=i*theta;
           A=[cos(theta), sin(theta), -1, 0; 
                 -sin(theta), cos(theta), 0, -1; 
                 -1, xi_1/(SRD+RDD), 0, 0,;
                 0, 0, -1, xi_2/(SRD+RDD)];
            b=[0;0;-xi_1*SRD/(SRD+RDD);-xi_2*SRD/(SRD+RDD);];
            x=(A\b);
            r0(i,1:2)=(T_2d(x(3:4), -alpha))';
            r0(i,3)=proj(i+1,2)*(SRD+x(4))/(RDD+SRD);
          
       end
   function [r2]=T_2d(r1,alpha)
   r2=[cos(-alpha) -sin(-alpha); sin(-alpha) cos(-alpha)]*r1;
   end
   end
   
   %%
   function [r0]= z_proj2r0(proj,theta,SRD,RDD)
   
   end

end