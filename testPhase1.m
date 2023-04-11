SRD = 1; % m, Source-Reference Distance
RDD = 1; % m, Reference-Detector (screen) Distance
theta_degree = 10; % clock-wise degree, what is the angle of camera rotation before each shot
NOS = [2:2:50]; % number of shots
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
plot(NOS(2:end),mean_average_deviation_percent(2:end))
legend(["Kalman" "Avg"])
title("Mean X deviation from true positions vs. NOS, 1000 runs per NOS");
xlabel("NOS");
ylabel("Mean deviation (%)")

