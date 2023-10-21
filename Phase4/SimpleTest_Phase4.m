%SOD = 15; ODD = 400;NOS=3; noise =0; delta_T = 1/68; %detector operates at 68FPS
% theta_degree=120*delta_T; %The tube turns 120degrees/sec
% theta=theta_degree*pi/180;
SOD = 1; ODD = 1;NOS=6; noise =0; delta_T = 1; theta=pi/6;
conditions=[noise,delta_T,NOS,theta,SOD,ODD];
r_expression=@(t) [10+4*t+t^2,20+5*t+t^2,6*t+30+t^2];
[xz_proj,real_positions] = generateTestPositions(r_expression, conditions);
predicted_values1=proj2r0_acc(xz_proj,theta,SOD,ODD,delta_T)
predicted_values2=proj2r0_acc_combination(xz_proj,theta,SOD,ODD,delta_T)
%add new particles
% r_expression1=@(t) [4*sin(t+0.1),4*cos(t-0.3),10*t];
% [xz_proj1,real_positions1] = generateTestPositions(r_expression1, conditions);
% plot(xz_proj(:,1),xz_proj(:,2)); title('projection plane'); hold on
% plot(xz_proj1(:,1),xz_proj1(:,2)); 
% figure
% plot3(real_positions(:,1),real_positions(:,2),real_positions(:,3))



