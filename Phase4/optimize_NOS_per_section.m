%% Find the best NOS_per_section based on the input
close all
clear

%% Independent Variables
rev=10; %revolutions of camera for the entire process
NL = 80; % lower bound
STP = 1;

% Input conditions
initial_positions = [0,0,0];
noise = 1e-3;
theta_degrees = 1.8;
camera_speed=1;%in Hz or revolution per second
weights = [0.5; 0.3; 0.2]; % Same order as distances

%% Variables that are not changed as frequent
SRD = 1; % m, Source-Reference Distance
RDD = 1; % m, Reference-Detector (screen) Distance
radius = 1;


%% Dependent Variables
NOS = rev*360/1.8
NH = 200; % higher bound
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
x = [NL:STP:NH]';
distances = [];

for i = 1:size(x,1)
    conditions(5) = x(i);
    %         distances(i) = computeMaxDistance(initial_positions,conditions,v,method,dataPiling,real_positions,xz_proj);
    [max_distance, sum_distance,methodTime] = computeMaxDistance(initial_positions,conditions,v,method,dataPiling,real_positions,xz_proj);
    distances = [distances;max_distance, sum_distance,methodTime];

end
weightedIndexes = distances*weights;
%% Plot
% Plot the data
figure('Position', [10 10 1200 800]) % set the figure size
figureNum = 4;

subplot(figureNum,1,1); % select the first subplot
[minVal, minIdx] = min(distances(:,1));
[maxVal, maxIdx] = max(distances(:,1));
plot(x,distances(:,1));
hold on;
plot((minIdx-1)*STP+x(1), minVal, 'ro'); % marks the min point with a red circle
plot((maxIdx-1)*STP+x(1), maxVal, 'bo'); % marks the max point with a blue circle
text((minIdx-1)*STP+x(1), minVal, sprintf('Min = %f, idx = %d', minVal,(minIdx-1)*STP+x(1)), 'VerticalAlignment','bottom');
text((maxIdx-1)*STP+x(1), maxVal, sprintf('Max = %f, idx = %d', maxVal,(maxIdx-1)*STP+x(1)), 'VerticalAlignment','top');
title('Max Distances');
hold off;

subplot(figureNum,1,2); % select the second subplot
[minVal, minIdx] = min(distances(:,2));
[maxVal, maxIdx] = max(distances(:,2));
plot(x,distances(:,2));
hold on;
plot((minIdx-1)*STP+x(1), minVal, 'ro'); % marks the min point with a red circle
plot((maxIdx-1)*STP+x(1), maxVal, 'bo'); % marks the max point with a blue circle
text((minIdx-1)*STP+x(1), minVal, sprintf('Min = %f, idx = %d', minVal,(minIdx-1)*STP+x(1)), 'VerticalAlignment','bottom');
text((maxIdx-1)*STP+x(1), maxVal, sprintf('Max = %f, idx = %d', maxVal,(maxIdx-1)*STP+x(1)), 'VerticalAlignment','top');
title('Sum Distances');
hold off;

subplot(figureNum,1,3); % select the second subplot
[minVal, minIdx] = min(distances(:,3));
[maxVal, maxIdx] = max(distances(:,3));
plot(x,distances(:,3));
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



% best_NOS_per_section = fminsearch(@(x)computeMaxDistance(x,initial_positions,conditions,v,method,dataPiling,real_positions,xz_proj), NOS_per_section_initial_guess);
% best_NOS_per_section


function [max_distance,sum_distance,methodTime] = computeMaxDistance(initial_positions,conditions,v,method,dataPiling,real_positions,xz_proj)
fprintf('currently on %d\n',conditions(5));

try
    %     conditions(5) = round(NOS_per_section);
    % Run your model
    tic
    estimated_positions = Phase4_trace_3d(initial_positions,conditions,v,method,dataPiling,xz_proj);
    methodTime = toc

    % Compute the distances
    distances = sqrt(sum((real_positions - estimated_positions).^2, 2));

    % Return the maximum distance
    max_distance = max(distances)
    sum_distance = sum(distances);

catch ME
    fprintf('An error occurred in Phase4_trace_3d: %s\n', ME.message);
    max_distance = inf;

end
end



% end
