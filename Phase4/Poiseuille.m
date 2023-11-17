% Poiseuille Flow Simulation with Analytical Solution in MATLAB

% Parameters
L = 1;          % Length of the pipe
W = 1;          % Width of the pipe
rho = 1;        % Density of the fluid
mu = 1;         % Viscosity of the fluid
dpdx = -1;      % Pressure gradient along the x-direction

% Grid
Nx = 50;        % Number of grid points along x
Ny = 50;        % Number of grid points along y
dx = L / Nx;    % Grid spacing along x
dy = W / Ny;    % Grid spacing along y

% Initialize velocity field
u = zeros(Nx+1, Ny+1);
v = zeros(Nx+1, Ny+1);


% Analytical solution for Poiseuille flow
R = W / 2;      % Radius of the pipe
y_analytical = linspace(-R, R, 100);
u_analytical = (1/(2*mu)) * (dpdx) * (R^2 - y_analytical.^2);

% Plot the analytical solution
figure;
plot(u_analytical, y_analytical, 'r-', 'LineWidth', 2);
xlabel('Velocity (u)');
ylabel('Radial Position (y)');
title('Analytical Solution for Poiseuille Flow');
grid on;


% Plot pathlines from numerical simulation
[x, y] = meshgrid(0:dx:L, 0:dy:W);
figure;
streamline(x, y, u, v, linspace(0, L, 20), linspace(0, W, 20));
xlabel('X-axis');
ylabel('Y-axis');
title('Poiseuille Flow - Numerical Simulation');
