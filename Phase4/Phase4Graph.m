function Phase4Graph(real_positions, estimated_positions, conditions,v,method,dataPiling,f2,animated,noise)
% close all

delta_T = conditions(2);
NOS = conditions(3);
theta_degree = conditions(4);
NOS_per_section = conditions(5);

real_positions = real_positions(1:height(estimated_positions),:);


figure(f2);
% % Get screen size
% screen = get(groot, 'ScreenSize');
% set(f2,'visible','on');
%
% % Set the new position [left bottom width height]
% % Increase the width of the figure by multiplying 'figure_outpos(3)' by a factor > 1
% % For example, let's increase the width by 50%:
% figure_outpos = get(gcf, 'OuterPosition');
% new_width = figure_outpos(3) * 1.5;  % Increase the width by 50%
% new_height = figure_outpos(4)*1.3;  % increase by 30%
%
% position = [(screen(3)-new_width)/2 (screen(4)-new_height)/2 new_width new_height];
%
% % Set the new position
% set(gcf, 'OuterPosition', position);

% Threshold (percent of the total distance traveled)
% threshold_percentage = 10;  % Set your own value here
% total_distance_travelled = sqrt(sum((real_Positions(1,:) - real_Positions(end,:)).^2));
% threshold = threshold_percentage/100 * total_distance_travelled; % threshold as a percentage of total distance travelled


% Calculate Euclidean distances between corresponding points
distances = sqrt(sum((real_positions - estimated_positions).^2, 2));

% Find indices of estimated points that are far from true points
% far_indices = find(distances > threshold);

% Find index of the minimum and maximum distance
[min_val, min_idx] = min(distances);
[max_val, max_idx] = max(distances);


% Define axes position [left, bottom, width, height] in normalized units
ax1 = axes(f2, 'Position', [0.1, 0.2, 0.6, 0.7]);  % move the plot down by 20%



% axes1 = axes(f2,)

if animated == true
    % Animation
    
    for k = 1:size(estimated_positions, 1)
        % Add the current estimated position to the plot

        
        plot3(ax1,real_positions(1:k, 1), real_positions(1:k, 2), real_positions(1:k, 3), 'r',"MarkerSize",10);
        hold on
        plot3(ax1,estimated_positions(1:k, 1), estimated_positions(1:k, 2), estimated_positions(1:k, 3), 'b', 'LineWidth', 2);

        % Add a text element showing the current NOS
        current_NOS_text = sprintf('Current NOS: %d', k);
        h = text(0.4, 0.8, current_NOS_text, 'Units', 'normalized');

        % Wait for a short time before the next frame
        pause(0.0001);

        % If it's not the last iteration, delete the previous line
        if k < size(estimated_positions, 1)
            children = get(gca, 'children');
            %             lines = findobj(gca, 'Type', 'line');

            if ~isempty(children)
                delete(children(1));
            end


            hold off
        end
    end
    delete(h);
    ax1.Position = [0.1, 0.2, 0.85, 0.7];


else
    %     Plotting and labeling
    ax1.Position = [0.1, 0.2, 0.85, 0.7];
    plot3(real_positions(:, 1), real_positions(:, 2), real_positions(:, 3), 'r',"MarkerSize",10);
    hold on;
    plot3(estimated_positions(:, 1), estimated_positions(:, 2), estimated_positions(:, 3), 'b', 'LineWidth', 2);

end

grid on;
xlabel(ax1,'X');
ylabel(ax1,'Y');
zlabel(ax1,'Z');
title(ax1,'Real Positions vs Estimated Positions in 3D');
legend(ax1,'Real Positions', 'Estimated Positions');
% %
% % Label the far points
% for i = 1:length(far_indices)
%     idx = far_indices(i);
%     text(estimatedPositions(idx, 1), estimatedPositions(idx, 2), estimatedPositions(idx, 3), sprintf('E%d', idx), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
% end

% Label the point with the minimum distance
text(real_positions(min_idx, 1), real_positions(min_idx, 2), real_positions(min_idx, 3), sprintf('Min (%d): (%0.2f, %0.2f, %0.2f)', min_idx, real_positions(min_idx, :)), 'Color', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(estimated_positions(min_idx, 1), estimated_positions(min_idx, 2), estimated_positions(min_idx, 3), sprintf('Min (%d): (%0.2f, %0.2f, %0.2f)', min_idx, estimated_positions(min_idx, :)), 'Color', 'b', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');

% Label the point with the maximum distance
text(real_positions(max_idx, 1), real_positions(max_idx, 2), real_positions(max_idx, 3), sprintf('Max (%d): (%0.2f, %0.2f, %0.2f)', max_idx, real_positions(max_idx, :)), 'Color', 'r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(estimated_positions(max_idx, 1), estimated_positions(max_idx, 2), estimated_positions(max_idx, 3), sprintf('Max (%d): (%0.2f, %0.2f, %0.2f)', max_idx, estimated_positions(max_idx, :)), 'Color', 'b', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');


annotation('textbox', [0.74, 0.7, 0.2, 0.1], 'String', 'Least Square Big Matrix Variable Acceleration', 'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'EdgeColor', 'none',FontWeight='bold');

% Adding the distances, threshold percentage and other information in a single annotation
annotation('textbox', [0.74, 0.2, 0.3, 0.2], 'String', ...
    sprintf(['Noise: %3f\n' ...
    'Min Distance: %0.2f\nMax Distance: %0.2f\n' ...
    'Number of Shots: %d\n' ...
    'Rotation: %d degrees\n' ...
    'Time Between Shots: %0.4f s\n' ...
    'NOS Per Section: %d\n' ...
    'Method: %s\n' ...
    'Data Piling: %s'], ...
    noise, min_val, max_val, NOS, theta_degree, delta_T,NOS_per_section,method,dataPiling), ...
    'FitBoxToText', 'on');

% Get the string representation of the velocity function
velocity_expression = func2str(v);

% Remove the '@(t)' part from the velocity expression string
velocity_expression = strrep(velocity_expression, '@(t)', '');

% Split the velocity expression into x, y, and z components
velocity_components = strsplit(velocity_expression, ',');

% Remove the brackets from each velocity component
velocity_components = strrep(velocity_components, '[', '');
velocity_components = strrep(velocity_components, ']', '');

% % Format the annotation string
% annotation_string = sprintf('True velocity:\nX: %s\nY: %s\nZ: %s\nt is time elapsed', velocity_components{1}, velocity_components{2}, velocity_components{3});
% 
% % Annotate function definition
% annotation('textbox', [0.74, 0.55, 0.26, 0.1], ...
%     'String', annotation_string, ...
%     'Units', 'normalized', ...
%     'HorizontalAlignment', 'left', ...
%     'VerticalAlignment', 'top', ...
%     'EdgeColor', 'none');



% Set the scale for each axis
x_scale = [-1, 1]; % Specify the desired minimum and maximum values for the x-axis
y_scale = [-1, 1]; % Specify the desired minimum and maximum values for the y-axis
% z_scale = [-1, 1]; % Specify the desired minimum and maximum values for the z-axis

% Set the scale of each axis
% xlim(x_scale);
% ylim(y_scale);
% zlim(z_scale);



% ... Your previous code ...

% Prompt for saving the figure
% saveImage(NOS, theta_degree, delta_T, f2);

end

function saveImage(NOS, theta_degree, delta_T, f2)
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
    saveas(f2, fullfile(folder_path, filename), 'png'); % You can change the file format if desired

    disp(['Figure saved as ', fullfile(folder_path, filename)]);
end
end