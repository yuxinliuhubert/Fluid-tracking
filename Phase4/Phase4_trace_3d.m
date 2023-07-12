%% Data Input
function positions_predicted = Phase4_trace_3d(varargin)
%N means the number of shots to use for one matrix; NOS/N must be a whole
%number
format default

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
addRequired(p,'method');
addRequired(p,'dataPiling');
addRequired(p,'xz_proj');

% parse inputs
parse(p,varargin{:});

% use the results
% initial_position_3d = p.Results.initial_position_3d;
conditions = p.Results.conditions;
% vel_expression = p.Results.vel_expression;
method = p.Results.method; % this will contain 'linear' or 'acceleration' depending on the input
dataPiling = p.Results.dataPiling;
xz_proj = p.Results.xz_proj;

[delta_T,NOS,theta_degree,N,SRD,RDD] = deal(conditions(2), ...
    conditions(3),conditions(4),conditions(5), ...
    conditions(6),conditions(7));

% method = 0; % 0 for least square, 1 for kalman

% v = @(t) [3*t+7;2;1*t+2.2]; %velocity function

% Process error
delta_P_X = 1e-3; % m
delta_P_VX = 1e-3; % m/sec



.....
    .....


% Observation error
delta_O_X= 2e-3; % m
delta_O_VX = 1.5e-4; % m
% delta_O_Y = delta_O_X;
% delta_O_VY = delta_O_VX;

%% Calculation part 1: projection generation
theta = theta_degree/180*pi;

%% Calculation part2: geting the measured values

proj_used_index=1;
% Preallocate positions_predicted
positions_predicted = zeros(NOS, 3);
NOS_per_section = N;
prev_NOS_section = NOS_per_section;

if strcmp(dataPiling,'serial')
    % no-overlap method
    while proj_used_index < NOS
        alpha=-theta*(proj_used_index-1);% alpha is for tracking the degree rotated from the 1st shot
        temp = generateEstimatedPositions(alpha, proj_used_index, NOS_per_section);
        positions_predicted(proj_used_index:proj_used_index+NOS_per_section-1, :)= temp;
%         count = count + size(temp,1);

        proj_used_index = proj_used_index+NOS_per_section; % change project_used_index accordingly
%             prev_NOS_section = NOS_per_section;
            if abs(NOS-proj_used_index) < N
                % Back track so that we have enough shots
                prev_proj_index = proj_used_index;
                proj_used_index = NOS - NOS_per_section +1; % Go back and redo
                alpha=-theta*(proj_used_index-1);% alpha is for tracking the degree rotated from the 1st shot
                last_positions = positions_predicted(proj_used_index: prev_proj_index-1, :);
                new_positions = generateEstimatedPositions(alpha, proj_used_index, NOS_per_section);
                combined_positions = [(new_positions(1:height(last_positions), :) + last_positions)/2; new_positions(height(last_positions)+1: end,:)];
                
                positions_predicted(proj_used_index:proj_used_index+NOS_per_section-1, :)= combined_positions;
                
                proj_used_index = proj_used_index + NOS_per_section+1;


%                 NOS_per_section = NOS - proj_used_index + 1; % adjust the NOS_per_section to include the remaining shots              
            end
%         if NOS - proj_used_index < N+6
%             prev_NOS_section = NOS_per_section;
%             NOS_per_section = NOS - proj_used_index + 1; % adjust the NOS_per_section to include the remaining shots
%         end
    end

elseif strcmp(dataPiling,'overlap')
    for i = 1:round(NOS-N)
        alpha=-theta*(proj_used_index-1);% alpha is for tracking the degree rotated from the 1st shot

        if proj_used_index == 1
            positions_predicted(proj_used_index:N,:) = generateEstimatedPositions(alpha,proj_used_index, NOS_per_section);
        else
            % take every N shots from every index, and take average of them
            last_positions = positions_predicted(proj_used_index: proj_used_index + prev_NOS_section-2, :);
            new_positions = generateEstimatedPositions(alpha,proj_used_index, NOS_per_section);

            % add the new positions with the old, and then take average
            new_positions = [(new_positions(1:height(last_positions), :) + last_positions)/2; new_positions(height(last_positions): end,:)];
            positions_predicted(proj_used_index:proj_used_index+size(new_positions,1)-1,:) = new_positions;
        end
        proj_used_index = proj_used_index + 1;
    end

end



    function positions_predicted = generateEstimatedPositions(alpha, proj_used_index, N)
    % Preallocate positions_predicted
    positions_predicted = zeros(N, 3);

    values_this_round=proj2r0_acc(xz_proj(proj_used_index:(proj_used_index+N-1),:),theta,SRD,RDD,delta_T);

    if strcmp(method,'acceleration')
        [x0, y0, z0 ,u ,v ,w ,a_x, a_y ,a_z]=deal(values_this_round(1),values_this_round(2),values_this_round(3),values_this_round(4),values_this_round(5),values_this_round(6),values_this_round(7),values_this_round(8),values_this_round(9));
        position_rotated=T([x0;y0;z0],alpha)';
        [x0, y0, z0]=deal(position_rotated(1), position_rotated(2) , position_rotated(3) );
        positions_predicted(1,:)=position_rotated;
        velocity_rotated=T([u;v;w],alpha)';
        [u ,v ,w]=deal(velocity_rotated(1), velocity_rotated(2) , velocity_rotated(3) );
        acc_rotated=T([a_x; a_y ;a_z],alpha)';
        [a_x, a_y ,a_z]=deal(acc_rotated(1), acc_rotated(2) , acc_rotated(3) );
        for j =2:N
            time=delta_T*(j-1);
            positions_predicted(j,:)=[x0+u*time+0.5*a_x*time^2, y0+v*time+0.5*a_y*time^2, z0+w*time+0.5*a_z*time^2];
        end
    elseif  strcmp(method,'linear')
        [x0, y0, z0 ,u ,v ,w]=deal(values_this_round(1),values_this_round(2),values_this_round(3),values_this_round(4),values_this_round(5),values_this_round(6));
        position_rotated=T([x0;y0;z0],alpha)';
        [x0, y0, z0]=deal(position_rotated(1), position_rotated(2) , position_rotated(3) );
        positions_predicted(1,:)=position_rotated;
        velocity_rotated=T([u;v;w],alpha)';
        [u ,v ,w]=deal(velocity_rotated(1), velocity_rotated(2) , velocity_rotated(3) );
        for j =2:N
            time=delta_T*(j-1);
            positions_predicted(j,:)=[x0+u*time, y0+v*time, z0+w*time];
        end
    end
end



end
