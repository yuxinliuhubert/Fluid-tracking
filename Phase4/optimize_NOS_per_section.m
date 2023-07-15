%% Find the best NOS_per_section based on the input
close all
clear

%% Independent Variables
rev=10; %revolutions of camera for the entire process

% Input conditions
initial_positions = [0,0,0];
noise = 1e-3;
theta_degrees = 1.8;
camera_speed=0.5;%in Hz or revolution per second
weights = [0.5; 0.3; 0.2]; % Same order as distances

%% Variables that are not changed as frequent
SRD = 1; % m, Source-Reference Distance
RDD = 1; % m, Reference-Detector (screen) Distance
radius = 1;


%% Dependent Variables
NOS = rev*360/1.8
% NH = 200; % higher bound
delta_T=camera_speed*theta_degrees/360;
shots_per_second = 1/delta_T;
% NOS=floor(360*rev/theta_degrees/NOS_per_section)*NOS_per_section;
v=@(t)[0.9*sin(t), 0.9*cos(t),1];

NOS_per_section_initial_guess = round(NOS/9);  % start the search from 120
conditions = [noise, delta_T, NOS,theta_degrees,NOS_per_section_initial_guess,SRD,RDD];
method = 'acceleration';
dataPiling = 'serial';

% Set up the conditions
[xz_proj, real_positions] = generateTestPositions(v,initial_positions,conditions);


%% Calculate the distances
% x = round(linspace(NL,NH,6))';
NL = 6;
NH = NOS/2;
STP = floor((NH-NL)/10);
if STP <= 1
    STP = 1;
end
fprintf('NL: %d, NH: %d, STP: %d\n',NL, NH,STP);

while STP >= 1
x = [NL:STP:NH]';
distances_info = [];

for i = 1:size(x,1)
    conditions(5) = x(i);
    %         distances(i) = computeMaxDistance(initial_positions,conditions,v,method,dataPiling,real_positions,xz_proj);
    [max_distance, sum_distance,methodTime] = computeMaxDistance(initial_positions,conditions,v,method,dataPiling,real_positions,xz_proj);
    distances_info = [distances_info;max_distance, sum_distance,methodTime];
    if methodTime > 20
        x = x(1:i);
        break;
    end

end
mean_distances_info = mean(distances_info,1);
weightedIndexes = (distances_info-mean_distances_info)./mean_distances_info *weights;

%% Plot
% Plot the data
figure('Position', [10 10 1200 800]) % set the figure size
figureNum = 4;

subplot(figureNum,1,1); % select the first subplot
[minVal, minIdx] = min(distances_info(:,1));
[maxVal, maxIdx] = max(distances_info(:,1));
plot(x,distances_info(:,1));
hold on;
plot((minIdx-1)*STP+x(1), minVal, 'ro'); % marks the min point with a red circle
plot((maxIdx-1)*STP+x(1), maxVal, 'bo'); % marks the max point with a blue circle
text((minIdx-1)*STP+x(1), minVal, sprintf('Min = %f, idx = %d', minVal,(minIdx-1)*STP+x(1)), 'VerticalAlignment','bottom');
text((maxIdx-1)*STP+x(1), maxVal, sprintf('Max = %f, idx = %d', maxVal,(maxIdx-1)*STP+x(1)), 'VerticalAlignment','top');
title('Max Distances');
hold off;

subplot(figureNum,1,2); % select the second subplot
[minVal, minIdx] = min(distances_info(:,2));
[maxVal, maxIdx] = max(distances_info(:,2));
plot(x,distances_info(:,2));
hold on;
plot((minIdx-1)*STP+x(1), minVal, 'ro'); % marks the min point with a red circle
plot((maxIdx-1)*STP+x(1), maxVal, 'bo'); % marks the max point with a blue circle
text((minIdx-1)*STP+x(1), minVal, sprintf('Min = %f, idx = %d', minVal,(minIdx-1)*STP+x(1)), 'VerticalAlignment','bottom');
text((maxIdx-1)*STP+x(1), maxVal, sprintf('Max = %f, idx = %d', maxVal,(maxIdx-1)*STP+x(1)), 'VerticalAlignment','top');
title('Sum Distances');
hold off;

subplot(figureNum,1,3); % select the second subplot
[minVal, minIdx] = min(distances_info(:,3));
[maxVal, maxIdx] = max(distances_info(:,3));
plot(x,distances_info(:,3));
hold on;
plot((minIdx-1)*STP+x(1), minVal, 'ro'); % marks the min point with a red circle
plot((maxIdx-1)*STP+x(1), maxVal, 'bo'); % marks the max point with a blue circle
text((minIdx-1)*STP+x(1), minVal, sprintf('Min = %f, idx = %d', minVal,(minIdx-1)*STP+x(1)), 'VerticalAlignment','bottom');
text((maxIdx-1)*STP+x(1), maxVal, sprintf('Max = %f, idx = %d', maxVal,(maxIdx-1)*STP+x(1)), 'VerticalAlignment','top');
title('Run Time');
hold off;

subplot(figureNum,1,4); % select the second subplot
[minVal, minIdx] = min(weightedIndexes);
[maxVal, maxIdx] = max(weightedIndexes);
plot(x,weightedIndexes);
hold on;
plot((minIdx-1)*STP+x(1), minVal, 'ro'); % marks the min point with a red circle
plot((maxIdx-1)*STP+x(1), maxVal, 'bo'); % marks the max point with a blue circle
text((minIdx-1)*STP+x(1), minVal, sprintf('Min = %f, idx = %d', minVal,(minIdx-1)*STP+x(1)), 'VerticalAlignment','bottom');
text((maxIdx-1)*STP+x(1), maxVal, sprintf('Max = %f, idx = %d', maxVal,(maxIdx-1)*STP+x(1)), 'VerticalAlignment','top');
title('Total Weighted');
hold off;

drawnow;

if minIdx == 1 || minIdx == length(x) || STP == 1
    STP = 0;
    break;
else
    NL = x(minIdx-1);
    NH = x(minIdx+1);
    STP = floor((NH-NL)/10);
    if STP <= 1
        STP = 1;

    end
end
fprintf('New search parameters: NL=%d, NH=%d, STP=%d\n',NL, NH,STP);

end

fprintf('Optimized NOS_per_Section: %d\n',(minIdx-1)*STP+x(1))

% best_NOS_per_section = fminsearch(@(x)computeMaxDistance(x,initial_positions,conditions,v,method,dataPiling,real_positions,xz_proj), NOS_per_section_initial_guess);
% best_NOS_per_section


function [max_distance,sum_distance,methodTime] = computeMaxDistance(initial_positions,conditions,v,method,dataPiling,real_positions,xz_proj)

try
    %     conditions(5) = round(NOS_per_section);
    % Run your model

    fprintf('currently on shot %d...', conditions(5))
    tic
    estimated_positions = Phase4_trace_3d(initial_positions,conditions,v,method,dataPiling,xz_proj);
    methodTime = toc;

    % Compute the distances
    distances = sqrt(sum((real_positions - estimated_positions).^2, 2));

    % Return the maximum distance
    max_distance = max(distances);
    sum_distance = sum(distances);
    fprintf('Max dist: %.3f m, Sum dist: %.3f m\n',max_distance,sum_distance);

catch ME
    fprintf('An error occurred in Phase4_trace_3d: %s\n', ME.message);
    max_distance = inf;

end
end



% end
