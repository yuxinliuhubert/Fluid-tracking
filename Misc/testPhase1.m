function [mean_kalman_deviation_percent, mean_average_deviation_percent] = testPhase1(theta_input, NOS_input, run_input, noise_input, position_x, position_y, position_z)
SRD = 1; % m, Source-Reference Distance
RDD = 1; % m, Reference-Detector (screen) Distance
theta_degree = theta_input; % clock-wise degree, what is the angle of camera rotation before each shot
NOS = [2:2:NOS_input]; % number of shots
delta_T = 1; % s, time between shots
noise= noise_input;

% positions
x = position_x;
y = position_y;
z = position_z;

true_position = [x y z]
true_mean_position = mean(true_position);
runs = run_input;
for j = 1:length(NOS)
fprintf("Calculating NOS = "+NOS(j)+"\n")
tic
for i = 1:runs
[kalman_position, avg_position] = Phase1_pt_3d(true_position, noise,delta_T, NOS(j), theta_degree, SRD, RDD);
% kalman_positions(i,j) = mean(kalman_position,"all");
% avg_positions(i,j) = mean(avg_position,"all");
kalman_positions(i,j) = kalman_position(1);
avg_positions(i,j) = avg_position(1);

end
toc
end

% kalman_deviations = kalman_positions - true_mean_position;
% average_deviations = avg_positions - true_mean_position;
kalman_deviations = abs(kalman_positions - true_position(1));
average_deviations = abs(avg_positions - true_position(1));

mean_kalman_deviation = mean(kalman_deviations);
mean_average_deviation = mean(average_deviations)

mean_kalman_deviation_percent = mean_kalman_deviation/true_position(1)*100
mean_average_deviation_percent = mean_average_deviation/true_position(1)*100
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

%% graph
% title(t, "Deviation from true position after " + runs+ " runs when", 'Fontweight','bold')
plot(NOS(2:end),mean_kalman_deviation_percent(2:end))
hold on
plot(NOS(1:end),mean_average_deviation_percent(1:end))
legend(["Kalman" "Avg"])
title("Mean X deviation from true positions vs. NOS, 1000 runs per NOS");
xlabel("NOS");
ylabel("Mean deviation (%)")

fprintf("\n running Phase 1 Cyrus code...")

end
