function Plotting(positions)
% Create the figure
figure;

NOS=size(positions); NOS=NOS(1);

% Plot the curve with gradually changing color
for i = 1:NOS-1
    % Compute the color based on the z-value
    col = [1-i/NOS, 0, i/NOS];
    x1=positions(i,1); y1=positions(i,2);z1=positions(i,3);
    x2=positions(i+1,1); y2=positions(i+1,2);z2=positions(i+1,3);
    % Draw a line segment with the computed color
    plot3([x1,x2], [y1, y2], [z1, z2], 'Color', col);
    hold on;
end

% Label the axes
xlabel('X');
ylabel('Y');
zlabel('Z');
legend('particle 1')
% Show the plot
hold off;
% plot3(positions(1:2,1),positions(1:2,2),positions(1:2,3),'Color',[1,0,0])
