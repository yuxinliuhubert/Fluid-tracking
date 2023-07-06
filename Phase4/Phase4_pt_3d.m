%% Data Input
function [real_positions,positions_predicted] = Phase4_pt_3d(initial_position_3d, conditions,vel_expression,method, dataPiling,varargin)
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


p = inputParser; % create parser object

% add required parameters
addRequired(p,'initial_position_3d');
addRequired(p,'conditions');
addRequired(p,'vel_expression');

% add optional parameter 'method' with default value 'linear'
addParameter(p,'method','linear',@(x) any(validatestring(x,{'linear','acceleration'})));
addParameter(p,'dataPiling','serial',@(x) any(validatestring(x,{'serial','overlap'})));

% parse inputs
parse(p,initial_position_3d,conditions,vel_expression,varargin{:});

% use the results
initial_position_3d = p.Results.initial_position_3d;
conditions = p.Results.conditions;
vel_expression = p.Results.vel_expression;
method = p.Results.method; % this will contain 'linear' or 'acceleration' depending on the input
dataPiling = p.Results.dataPiling;

[noise, delta_T,NOS,theta_degree,N] = deal(conditions(1),conditions(2),conditions(3),conditions(4),conditions(5));

% method = 0; % 0 for least square, 1 for kalman

% v = @(t) [3*t+7;2;1*t+2.2]; %velocity function

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
    r0_k=r0_0+integral(vel_expression,0,delta_T*k,'ArrayValued', true); % true location in the original frame of reference
    real_positions = [real_positions; r0_k];
    r0_k_rotated=T(r0_k',theta*k)';
    M_p = (SRD+RDD)/(SRD+r0_k_rotated(2)); % magnification of particle
    x_proj=[x_proj;M_p*r0_k_rotated(1)+randn*noise/2];
    z_proj=[z_proj;M_p*r0_k_rotated+randn*noise/2];
end
% x_proj and z_proj are NOS by 1 vectors with all the x-coordinates projection on the screen
xz_proj=[x_proj, z_proj];

proj_used_index=1;
positions_predicted=[];
%% Calculation part2: geting the measured values

NOS_per_section = N;
prev_NOS_section = NOS_per_section;


if strcmp(dataPiling,'serial')
    % % no-overlap method
    for i = 1:round(NOS/N)-1

        alpha=-theta*(proj_used_index-1);% alpha is for tracking the degree rotated from the 1st shot
        positions_predicted = [positions_predicted;generateEstimatedPositions(alpha, proj_used_index, NOS_per_section)]
        % when mod(NOS, N) is not 0, this code ensures that at least 5
        % shots are used to generate positionsproj_used_index
        if NOS - proj_used_index < 2*N % insufficient equation trigger
            prev_NOS_section = NOS_per_section;
            NOS_per_section = NOS - proj_used_index; % adjust the NOS_per_section to include the remaining shots
        end
        proj_used_index = proj_used_index+NOS_per_section; % change project_used_index accordingly
    end

elseif strcmp(dataPiling,'overlap')
    for i = 1:round(NOS-N)
        alpha=-theta*(proj_used_index-1);% alpha is for tracking the degree rotated from the 1st shot

        if proj_used_index == 1
            positions_predicted = [positions_predicted; generateEstimatedPositions(alpha,proj_used_index, NOS_per_section)];
        else
            % take every N shots from every index, and take average of them
            last_positions = positions_predicted(proj_used_index: proj_used_index + prev_NOS_section-2, :);
            new_positions = generateEstimatedPositions(alpha,proj_used_index, NOS_per_section);

            % add the new positions with the old, and then take average
            new_positions = [(new_positions(1:height(last_positions), :) + last_positions)/2; new_positions(height(last_positions): end,:)];
            positions_predicted = [positions_predicted(1:proj_used_index-1,:);new_positions];
        end
        proj_used_index = proj_used_index + 1;
    end

end

%% function for rotating coordinate axes by a given angle alpha
    function [r2]=T(r1,alpha) %r1=[x1;y1;z1]and r2=[x2;y2;z2] describe the same physical pt but in 2 coordinate systems with
        % theta degress rotation CCW
        r2=[cos(-alpha) -sin(-alpha) 0; sin(-alpha) cos(-alpha) 0;0 0 1]*r1;
        % note that we are rotating the coordinate axes, not the location vector of the particle
    end



    function positions_predicted = generateEstimatedPositions(alpha, proj_used_index, N)
        positions_predicted = [];
        values_this_round=proj2r0_acc(xz_proj(proj_used_index:(proj_used_index+N-1),:),theta,SRD,RDD,delta_T);

        if strcmp(method,'acceleration')
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
        elseif  strcmp(method,'linear')
            [x0, y0, z0 ,u ,v ,w]=deal(values_this_round(1),values_this_round(2),values_this_round(3),values_this_round(4),values_this_round(5),values_this_round(6));
            position_rotated=T([x0;y0;z0],alpha)';
            [x0, y0, z0]=deal(position_rotated(1), position_rotated(2) , position_rotated(3) );
            positions_predicted=[positions_predicted;position_rotated];
            velocity_rotated=T([u;v;w],alpha)';
            [u ,v ,w]=deal(velocity_rotated(1), velocity_rotated(2) , velocity_rotated(3) );
            for j =1:N-1
                time=delta_T*j;
                positions_predicted=[positions_predicted;x0+u*time, y0+v*time, z0+w*time];
            end

        end

    end


end
