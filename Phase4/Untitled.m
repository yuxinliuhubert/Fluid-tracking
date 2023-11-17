% Some random data
x = linspace(0, 10, 200);
y = sin(x);

% A variable to represent color gradient
z = x;  % Color gradient based on x values

% Create a scatter plot with a color gradient
scatter(x, y, 50, z, 'filled');

% % Add a colorbar to show the color scale
% colorbar;
% 
% % Apply a colormap
% colormap(jet);
