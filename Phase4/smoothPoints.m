function filtered_measurements = smoothPoints(estimated_positions1,method,framesize)

if strcmp(method,'avg')
% Define the window size for the moving average
windowSize = framesize;

% Create a moving average filter
b = (1/windowSize)*ones(1, windowSize);
a = 1;

% Apply the moving average filter to each dimension
filtered_measurements = estimated_positions1;
filtered_measurements(:,1) = filter(b, a, estimated_positions1(:,1));
filtered_measurements(:,2) = filter(b, a, estimated_positions1(:,2));
filtered_measurements(:,3) = filter(b, a, estimated_positions1(:,3));

elseif strcmp(method,'sg')
    % Your 3D measurement data as a Nx3 matrix
measurements = estimated_positions1; 

% Define the polynomial order for the Savitzky-Golay filter
polyorder = 2;

% Apply the Savitzky-Golay filter to each dimension
filtered_measurements = measurements;
filtered_measurements(:,1) = sgolayfilt(measurements(:,1), polyorder, framesize);
filtered_measurements(:,2) = sgolayfilt(measurements(:,2), polyorder, framesize);
filtered_measurements(:,3) = sgolayfilt(measurements(:,3), polyorder, framesize);

elseif strcmp(method,'cb')
    % Your 3D measurement data as a Nx3 matrix
measurements = estimated_positions1;

% Define the smoothing parameter for csaps
% The value should be between 0 and 1. Higher values result in less smoothing.
p = 0.5; 

% Generate the array for the x-values assuming equally spaced points
x = 1:size(measurements, 1);

% Apply the csaps smoothing to each dimension
filtered_measurements = measurements;
filtered_measurements(:,1) = csaps(x, measurements(:,1), p, x);
filtered_measurements(:,2) = csaps(x, measurements(:,2), p, x);
filtered_measurements(:,3) = csaps(x, measurements(:,3), p, x);

end



end