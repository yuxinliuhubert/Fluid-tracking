%% Data Input
function results = Phase3_pt_3d_success(initial_position_3d, noise, delta_T, NOS,theta_degree)
% clear figure
format short
close all
SRD = 1; % m, Source-Reference Distance
RDD = 1; % m, Reference-Detector (screen) Distance
% % theta_degree = 10; % clock-wise degree, what is the angle of camera rotation before each shot
% % NOS = 100; % number of shots
% % delta_T = time between shots
% % initial_position_3d = [x0;y0;z0], a column vector
% noise=5e-3; %the 68% chance deviation of the projection measurement from the perfect measurement 

method = 0; % 0 for least square, 1 for kalman

% v = @(t) [3*t^2+7;2;1*t+2.2*t^2]; %velocity function
v = @(t) [5;1;1]; %velocity function

poly_highest_degree = 2;
% v = @(t) [3;2;1]; %velocity function

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
x_proj=M_p*r0_0(1)+randn*noise/2;
z_proj=M_p*r0_0(3)+randn*noise/2;
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


results = proj2r0(xz_proj,theta,SRD,RDD,delta_T);

init_pos_est = results(1:3);
init_vel_est = results(4:6);
init_acc_est = results(7:9);
estimatedPositions = init_pos_est;

for k = 1:NOS-1
    current_time = delta_T*k;
    estimatedPositions = [estimatedPositions; init_pos_est+ init_vel_est*current_time + 0.5*init_acc_est*current_time^2];

end

% plot the true positions in red ('r') and estimated positions in blue ('b')
% Create a figure
f1 = figure;


% Get screen size
screen = get(groot, 'ScreenSize');

% Set the new position [left bottom width height]
% Increase the width of the figure by multiplying 'figure_outpos(3)' by a factor > 1
% For example, let's increase the width by 50%:
figure_outpos = get(gcf, 'OuterPosition');
new_width = figure_outpos(3) * 1.5;  % Increase the width by 50%
new_height = figure_outpos(4)*1.3;  % increase by 30%

position = [(screen(3)-new_width)/2 (screen(4)-new_height)/2 new_width new_height];

% Set the new position
set(gcf, 'OuterPosition', position);

% Threshold (percent of the total distance traveled)
% threshold_percentage = 10;  % Set your own value here
% total_distance_travelled = sqrt(sum((real_Positions(1,:) - real_Positions(end,:)).^2));
% threshold = threshold_percentage/100 * total_distance_travelled; % threshold as a percentage of total distance travelled


% Calculate Euclidean distances between corresponding points
distances = sqrt(sum((real_Positions - estimatedPositions).^2, 2));

% Find indices of estimated points that are far from true points
% far_indices = find(distances > threshold);

% Find index of the minimum and maximum distance
[min_val, min_idx] = min(distances);
[max_val, max_idx] = max(distances);


% Define axes position [left, bottom, width, height] in normalized units
% ax = axes(f1, 'Position', [0.1, 0.2, 0.8, 0.6]);  % move the plot down by 20%

% Plotting and labeling
plot3(real_Positions(:, 1), real_Positions(:, 2), real_Positions(:, 3), 'r', 'LineWidth', 2);
hold on;
plot3(estimatedPositions(:, 1), estimatedPositions(:, 2), estimatedPositions(:, 3), 'b', 'LineWidth', 2);
% 
% % Label the far points
% for i = 1:length(far_indices)
%     idx = far_indices(i);
%     text(estimatedPositions(idx, 1), estimatedPositions(idx, 2), estimatedPositions(idx, 3), sprintf('E%d', idx), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
% end

% Label the point with the minimum distance
text(real_Positions(min_idx, 1), real_Positions(min_idx, 2), real_Positions(min_idx, 3), sprintf('Min (%d): (%0.2f, %0.2f, %0.2f)', min_idx, real_Positions(min_idx, :)), 'Color', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(estimatedPositions(min_idx, 1), estimatedPositions(min_idx, 2), estimatedPositions(min_idx, 3), sprintf('Min (%d): (%0.2f, %0.2f, %0.2f)', min_idx, estimatedPositions(min_idx, :)), 'Color', 'b', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');

% Label the point with the maximum distance
text(real_Positions(max_idx, 1), real_Positions(max_idx, 2), real_Positions(max_idx, 3), sprintf('Max (%d): (%0.2f, %0.2f, %0.2f)', max_idx, real_Positions(max_idx, :)), 'Color', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(estimatedPositions(max_idx, 1), estimatedPositions(max_idx, 2), estimatedPositions(max_idx, 3), sprintf('Max (%d): (%0.2f, %0.2f, %0.2f)', max_idx, estimatedPositions(max_idx, :)), 'Color', 'b', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');


annotation('textbox', [0.74, 0.7, 0.2, 0.1], 'String', 'Least Square Big Matrix Constant Acceleration', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'EdgeColor', 'none',FontWeight='bold');

% Adding the distances, threshold percentage and other information in a single annotation
annotation('textbox', [0.74, 0.2, 0.3, 0.2], 'String', ...
    sprintf(['Min Distance: %0.2f\nMax Distance: %0.2f\n' ...
             'Number of Shots: %d\n' ...
             'Rotation: %d degrees\n' ...
             'Time Between Shots: %0.2f s'], ...
             min_val, max_val, NOS, theta_degree, delta_T), ...
    'FitBoxToText', 'on');

% Get the string representation of the velocity function
velocity_expression = func2str(v);

% Remove the '@(t)' part from the velocity expression string
velocity_expression = strrep(velocity_expression, '@(t)', '');

% Split the velocity expression into x, y, and z components
velocity_components = strsplit(velocity_expression, ';');

% Remove the brackets from each velocity component
velocity_components = strrep(velocity_components, '[', '');
velocity_components = strrep(velocity_components, ']', '');

% Format the annotation string
annotation_string = sprintf('True velocity:\nX: %s\nY: %s\nZ: %s\nt is time elapsed', velocity_components{1}, velocity_components{2}, velocity_components{3});

% Annotate function definition
annotation('textbox', [0.74, 0.55, 0.26, 0.1], ...
    'String', annotation_string, ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', ...
    'EdgeColor', 'none');


grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Real Positions vs Estimated Positions in 3D');
legend('Real Positions', 'Estimated Positions');

% ... Your previous code ...

% Prompt for saving the figure
prompt = 'Do you want to save the figure? (y/n): ';
user_input = input(prompt, 's');

if strcmpi(user_input, 'y')
    % Define the format for the filename
    filename_format = 'NOS%d_R%d_T%.2f_D%d.png';
    
    % Replace the placeholders with appropriate values
    filename = sprintf(filename_format, NOS, theta_degree, delta_T, poly_highest_degree);
    
    % Specify the folder path
    folder_path = 'Generated Images/';
    
    % Save the figure in the specified folder
    saveas(f1, fullfile(folder_path, filename), 'png'); % You can change the file format if desired
    
    disp(['Figure saved as ', fullfile(folder_path, filename)]);
end




 %% function for rotating coordinate axes by a given angle alpha 
   function [r2]=T(r1,alpha) %r1=[x1;y1;z1]and r2=[x2;y2;z2] describe the same physical pt but in 2 coordinate systems with
   % theta degress rotation CCW
   r2=[cos(-alpha) -sin(-alpha) 0; sin(-alpha) cos(-alpha) 0;0 0 1]*r1;
   % note that we are rotating the coordinate axes, not the location vector of the particle
   end
  
  
%row_number_A=row_number_A+NOS ;%which is different from phase 1. We get a new equation each shot for the change of z coordinate due to velocity

function [r0]= proj2r0(proj,theta,SRD,RDD,delta_T) %takes in any N (even #) by 2 of projection coordinates and the angle 
   %theta btw the shots and outputs the predicted [x0, y0, z0] sets in coordinates when theta=0 (original
   %position before camera machine starts rotating
    NOS=height(proj); SDD=(SRD+RDD);
   row_number_A=round(2*NOS+ 2*(NOS-1)  ,  0  ); %2*NOS equations are from magnifications, and 
   % 2*(NOS-1) equations are from transformation
   col_number_A=round(1+2*NOS, 0); %round function is for avoiding precision error
   A=zeros(   row_number_A,  col_number_A); 
   b=zeros(row_number_A,1);
       for j = 1:(NOS)  %This for loop is for constructing the equations arising from magnification alone
           %for each increased number of shots there are 2 new variables
           %introduced and 2 equations
           xi_j=proj(j,1);   zi_j=proj(j,2); 
           A(2*j-1,1)=1; A( 2*j-1,2*j+1)=-zi_j/SDD; b(2*j-1)=zi_j*SRD/SDD;%for z0 magnification eq
           A(2*j,2*j)=-1;A(2*j,2*j+1)=xi_j/SDD; b(2*j)=-xi_j*SRD/SDD;%for x_0 magnification eq
           %the following 2 rows of codes are for equations transformation extracted from transformation
       end
       
       x=2*NOS+1; %x is for tracking the index of unfilled rows of the big matrix A
      for k = 2:(NOS)
          A( x:x+1, 2*k : 2*k+1 )=[-1 0; 0 -1];
          A(x:x+1,2:3)=[cos(theta*(k-1)) sin(theta*(k-1)); -sin(theta*(k-1)) cos(theta*(k-1))];
          
          x=x+2;
      end
       
%% Now, we expand the number of columns to incorporate new variables: u, v , w, a_x, a_y, a_z
A=[A,zeros(height(A),6)];
new_col_num=size(A,2);%new_col_num is the number of columns after adding the new variables 
%u, v , w, a_x, a_y, a_z are at the last 6 columns

%% The following is for coefficients related to V and a to magnification equations 
     for j = 1:(NOS)  
           %w is at new_col_num (the last column)
           A( 2*j-1,new_col_num-3)=delta_T*(j-1);
           A( 2*j-1,new_col_num)=0.5*(delta_T*(j-1))^2;%for w term 
           
     end
     % The following is for adding equations related to u, v, a_x, a_y to
     % transformation equations
     IoR=2*NOS+1; %IoR is for tracking the index of unfilled rows of the big matrix A
     
     for k=2:NOS  
         
         NoT=k-1; %NoT stands for the number of intervals, each with delta_T and theta, between the shot being
         % transformed and the first shot
              A(IoR, new_col_num-5)=cos(theta*NoT)*delta_T*(NoT);
              A(IoR+1, new_col_num-5)=-sin(theta*NoT)*delta_T*(NoT);%for u
              A(IoR, new_col_num-4)=sin(theta*NoT)*delta_T*(NoT);
              A(IoR+1, new_col_num-4)=cos(theta*NoT)*delta_T*(NoT);%for v
              A(IoR, new_col_num-2)=0.5*cos(theta*NoT)*(delta_T*NoT)^2;
              A(IoR+1, new_col_num-2)=-sin(theta*NoT)*0.5*(delta_T*NoT)^2;%for ax
              A(IoR, new_col_num-1)=0.5*sin(theta*NoT)*(delta_T*NoT)^2;
              A(IoR+1, new_col_num-1)=0.5*cos(theta*NoT)*(delta_T*NoT)^2;%for ay

             
            
         
          IoR=IoR+2;
     end

     x=(A\b);
r0=[x(2),x(3),x(1),x(new_col_num-5),x(new_col_num-4),x(new_col_num-3),x(new_col_num-2),x(new_col_num-1),x(new_col_num)];
   end

   

end