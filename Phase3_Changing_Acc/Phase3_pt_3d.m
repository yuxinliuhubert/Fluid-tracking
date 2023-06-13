%% Data Input
function results = Phase3_pt_3d(initial_position_3d, noise, delta_T, NOS,theta_degree)
format short
SRD = 1; % m, Source-Reference Distance
RDD = 1; % m, Reference-Detector (screen) Distance
% % theta_degree = 10; % clock-wise degree, what is the angle of camera rotation before each shot
% % NOS = 100; % number of shots
% % delta_T = time between shots
% % initial_position_3d = [x0;y0;z0], a column vector
% noise=5e-3; %the 68% chance deviation of the projection measurement from the perfect measurement 

method = 0; % 0 for least square, 1 for kalman

v = @(t) [0.1;t+0.1;2*t+0.2]; %velocity function

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
real_Positions = initial_position_3d;
M_p = (SRD+RDD)/(SRD+r0_0(2));
x_proj=[M_p*r0_0(1)+randn*noise/2];
z_proj=[M_p*r0_0(3)+randn*noise/2];
% % Second shot: with theta degrees rotation each time
for k = 1:NOS-1
    r0_k=r0_0+integral(v,0,delta_T*k,'ArrayValued', true); % true location in the original frame of reference
    real_Positions = [real_Positions; r0_0 + integral(v,0,delta_T*k,'ArrayValued', true)'];
    r_now=T(r0_k,theta*k);
    M_p = (SRD+RDD)/(SRD+r_now(2)); % magnification of particle
    x_proj=[x_proj;M_p*r_now(1)+randn*noise/2];
    z_proj=[z_proj;M_p*r_now(3)+randn*noise/2];
end
% x_proj and z_proj are NOS by 1 vectors with all the x-coordinates projection on the screen
xz_proj=[x_proj, z_proj];



%% Calculation part2: geting the measured values 

if method == 0
proj2r0(xz_proj,theta,SRD,RDD,delta_T)




else
xyz_measured=proj2r0(xz_proj,theta,SRD,RDD);

% initialize the states
xyz_initial = xyz_measured(1,:);
xyz_predicted(1,:) = xyz_initial;
velocity_predicted(1,:) = [0, 0, 0]; % x, y, z

% a state involves position and velocity. x_ki or x_kp will be a 2-by-1
% vector
% initial states
t= 0;
v_initial = v(0);
x_ki = [xyz_measured(1,1); v_initial(1)];
y_ki = [xyz_measured(1,2); v_initial(2)];
z_ki = [xyz_measured(1,3); v_initial(3)];
    
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

outputPosition = [x_k_prev,y_k_prev,z_k_prev];
t = t + delta_T;
if (NOS > 2)
% iteration begins
for i = 2:height(xyz_measured) 

    a_current = a(t);
    
    % predicted states
    x_kp = A * x_k_prev + B * a_current(1);
    y_kp = A * y_k_prev + B * a_current(2);
    z_kp = A * z_k_prev + B * a_current(3);

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
    outputPosition = [outputPosition; x_k_prev,y_k_prev,z_k_prev];
    t = t+delta_T;

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


end



%% Part 3: System output
% kalman_position_velocity = [x_k y_k z_k];
% r0_avg=mean(xyz_measured);

kalman_positions = outputPosition;

% Assign each variable to a column in the table
results = table(kalman_positions(:,1), real_Positions(:,1), kalman_positions(:,2), real_Positions(:,2), kalman_positions(:,3), real_Positions(:,3),...
    'VariableNames', {'kalman_x', 'real_x', 'kalman_y', 'real_y', 'kalman_z', 'real_z'});

% Pad xyz_measured if it's shorter
length_xyz = size(xyz_measured, 1);
length_kalman = size(kalman_positions, 1);
figure;

% Create x-axis for measured data
x_measured = 1:2:length_kalman;

% Plot Kalman x vs Real x vs Measured x
subplot(3,1,1);
plot(results.kalman_x, 'b', 'LineWidth', 1.5); % plot in blue
hold on;
plot(results.real_x, 'r', 'LineWidth', 1.5); % plot in red
plot(x_measured, xyz_measured(:,1), 'g', 'LineWidth', 1.5); % plot in green
hold off;
legend('Kalman x', 'Real x', 'Measured x');
ylabel('x-coordinate');
title('Comparison of Kalman x, Real x, and Measured x');

% Plot Kalman y vs Real y vs Measured y
subplot(3,1,2);
plot(results.kalman_y, 'b', 'LineWidth', 1.5); % plot in blue
hold on;
plot(results.real_y, 'r', 'LineWidth', 1.5); % plot in red
plot(x_measured, xyz_measured(:,2), 'g', 'LineWidth', 1.5); % plot in green
hold off;
legend('Kalman y', 'Real y', 'Measured y');
ylabel('y-coordinate');
title('Comparison of Kalman y, Real y, and Measured y');

% Plot Kalman z vs Real z vs Measured z
subplot(3,1,3);
plot(results.kalman_z, 'b', 'LineWidth', 1.5); % plot in blue
hold on;
plot(results.real_z, 'r', 'LineWidth', 1.5); % plot in red
plot(x_measured, xyz_measured(:,3), 'g', 'LineWidth', 1.5); % plot in green
hold off;
legend('Kalman z', 'Real z', 'Measured z');
xlabel('Time');
ylabel('z-coordinate');
title('Comparison of Kalman z, Real z, and Measured z');

% Improve the spacing between the plots
subplot_tool = subplot(3,1,3);
subplot_tool.Position(2) = subplot_tool.Position(2) - 0.05;
subplot_tool = subplot(3,1,2);
subplot_tool.Position(2) = subplot_tool.Position(2) - 0.025;




 %% function for rotating coordinate axes by a given angle alpha 
   function [r2]=T(r1,alpha) %r1=[x1;y1;z1]and r2=[x2;y2;z2] describe the same physical pt but in 2 coordinate systems with
   % theta degress rotation CCW
   r2=[cos(-alpha) -sin(-alpha) 0; sin(-alpha) cos(-alpha) 0;0 0 1]*r1;
   % note that we are rotating the coordinate axes, not the location vector of the particle
   end
  
   %% 
%row_number_A=row_number_A+NOS ;%which is different from phase 1. We get a new equation each shot for the change of z coordinate due to velocity

function [r0]= proj2r0(proj,theta,SRD,RDD,delta_T) %takes in any N (even #) by 2 of projection coordinates and the angle 
   %theta btw the shots and outputs the predicted [x0, y0, z0] sets in coordinates when theta=0 (original
   %position before camera machine starts rotating
   NOS=height(proj); SDD=(SRD+RDD);
   row_number_A=round(2*NOS+ 2*(factorial(NOS)/(factorial(NOS-2)*2))  ,  0  ); 
   col_number_A=round(1+2*NOS, 0); %round function is for avoiding precision error
   A=zeros(   row_number_A,  col_number_A); 
   b=zeros(height(A),1);
       for j = 1:(NOS)  %This for loop is for constructing the equations arising from magnification alone
           %for each increased number of shots there are 2 new variables
           %introduced and 2 equations
           xi_j=proj(j,1);   zi_j=proj(j,2); 
           A(2*j-1,1)=1; A( 2*j-1,2*j+1)=-zi_j/SDD; b(2*j-1)=zi_j*SRD/SDD;%for z0 magnification eq
           A(2*j,2*j)=-1;A(2*j,2*j+1)=xi_j/SDD; b(2*j)=-xi_j*SRD/SDD;%for x_0 magnification eq
           %the following 2 rows of codes are for equations transformation extracted from transformation
       end
       %Now, 1 to NOS*2 rows are filled
       IoR=2*NOS+1; %NoE, standing for Index of Rows, is for tracking the index of unfilled rows 
       %of the big matrix A
      for k = 2:(NOS)
          A( IoR:IoR+2*(k-1)-1, 2*k : 2*k+1 )=repmat([-1 0; 0 -1],k-1,1);
          
          for L=1:k-1
              A(IoR+2*(L-1) : IoR+2*(L-1)+1 , 2*(k-1)-2*(L-1):2*(k-1)-2*(L-1)+1)=[cos(theta*L) sin(theta*L); -sin(theta*L) cos(theta*L)];
          end
          IoR=IoR+2*(k-1);
      end
       
%Now, we expand the number of columns to incorporate new variables: u, v , w, a_x, a_y, a_z 
A=[A,zeros(height(A),6)];
new_col_num=length(A); %new_col_num is the number of columns after adding the new variables related to v and a 
%u, v , w, a_x, a_y, a_z are at the last 6 columns
%% The following is for coefficients related to V and a to magnification equations 
     for j = 1:(NOS-1)  
           %w is at new_col_num-3; a_z is at new_col_num (the last column)
           A( 2*j-1,new_col_num-3)=delta_T*(j-1); %for w term 
           A( 2*j-1,new_col_num)=0.5*(delta_T*(j-1))^2; %for a_z term 
     end
     % The following is for adding equations related to u, v, a_x, a_y to
     % transformation equations
     IoR=2*NOS+1; %IoR is for tracking the index of unfilled rows of the big matrix A
     for k=2:NOS  
         
          for L=1:k-1
              A(IoR+2*(L-1), new_col_num-5)=cos(theta*L)*delta_T*(L-1);
              A(IoR+2*(L-1)+1, new_col_num-5)=-sin(theta*L)*delta_T*(L-1);%for u
              A(IoR+2*(L-1), new_col_num-4)=sin(theta*L)*delta_T*(L-1);
              A(IoR+2*(L-1)+1, new_col_num-4)=cos(theta*L)*delta_T*(L-1);%for v
              A(IoR+2*(L-1), new_col_num-2)=0.5*cos(theta*L)*(delta_T*(L-1))^2;
              A(IoR+2*(L-1)+1, new_col_num-2)=-sin(theta*L)*0.5*(delta_T*(L-1))^2;%for ax
              A(IoR+2*(L-1), new_col_num-1)=0.5*sin(theta*L)*(delta_T*(L-1))^2;
              A(IoR+2*(L-1)+1, new_col_num-1)=0.5*cos(theta*L)*(delta_T*(L-1))^2;%for ay

          end
          IoR=IoR+2*(k-1);
     end

     x=(A\b);
     r0=[x(2),x(3),x(1),x(new_col_num-5),x(new_col_num-4),x(new_col_num-3),x(new_col_num-2),x(new_col_num-1),x(new_col_num)];
   end

   

end