SRD = 1; % m, Source-Reference Distance
RDD = 1; % m, Reference-Detector (screen) Distance
theta_degree = 10; % clock-wise degree, what is the angle of camera rotation before each shot
NOS = [10,100,200,500,800,1000,2000,5000]; % number of shots
delta_T = 1; % s, time between shots
noise= 5e-3;

% positions
x = 0.1;
y = 0.2;
z = 0.1;

true_position = [x y z]
true_mean_position = mean(true_position);
runs = 1000;
for j = 1:length(NOS)
fprintf("Calculating NOS = "+NOS(j)+"\n")
tic
for i = 1:runs
[kalman_position, avg_position] = Phase1_pt_3d(true_position, noise,delta_T, NOS(j), theta_degree, SRD, RDD);
% kalman_positions(i,j) = mean(kalman_position,"all");
% avg_positions(i,j) = mean(avg_position,"all");
kalman_positions(i,j) = kalman_position(3);
avg_positions(i,j) = avg_position(3);

end
toc
end

% kalman_deviations = kalman_positions - true_mean_position;
% average_deviations = avg_positions - true_mean_position;
kalman_deviations = kalman_positions - true_position(3);
average_deviations = avg_positions - true_position(3);

% t = tiledlayout(3,2);
% title(t, "Deviation from true position after " + runs+ " runs", 'Fontweight','bold')
% % x coordinates deviation compare
% nexttile
% histogram(kalman_deviations(:,1))
% title('x, kalman')
% 
% nexttile
% histogram(average_deviations(:,1))
% title('x, avg')
% 
% % y coordinates deviation compare
% nexttile
% histogram(kalman_deviations(:,2))
% title('y, kalman')
% 
% nexttile
% histogram(average_deviations(:,2))
% title('y, avg')
% 
% 
% % z coordinates deviation compare
% nexttile
% histogram(kalman_deviations(:,3))
% title('z, kalman')
% 
% nexttile
% histogram(average_deviations(:,3))
% title('z, avg')

t = tiledlayout(6,4);
% title(t, "Deviation from true position after " + runs+ " runs when", 'Fontweight','bold')
title(t, "Deviation from true z position after " + runs+ " runs when", 'Fontweight','bold')
for k = 1:length(NOS) 
    nexttile
    histogram(kalman_deviations(:,k))
    title("NOS = "+NOS(k) +", kalman")
    nexttile
    histogram(average_deviations(:,k))
    title("NOS = "+NOS(k) +", avg")
end



