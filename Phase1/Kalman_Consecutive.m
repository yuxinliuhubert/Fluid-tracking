function [percent_error]=Kalman_Consecutive(true_position_3d,noise,NOS, theta_degree)
delta_T=1; SRD=1;RDD=1;
[kalman_position, avg_position] = Phase1_pt_3d(true_position_3d, noise,delta_T, NOS, theta_degree, SRD, RDD);
percent_error=abs((avg_position(1)-true_position_3d(1))/true_position_3d(1))*100;