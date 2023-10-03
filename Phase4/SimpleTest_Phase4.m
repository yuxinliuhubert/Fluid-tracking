SOD = 15; ODD = 400;NOS=180; noise =0.1; delta_T = 1/68; %detector operates at 68FPS

theta_degree=120*delta_T; %The tube turns 120degrees/sec

theta=theta_degree*pi/180;
conditions=[noise,delta_T,NOS,theta_degree,SOD,ODD];
r_expression=@(t) [5*sin(t),5*cos(t),10*t];
[xz_proj,real_positions] = generateTestPositions(r_expression, conditions);
predicted_values=proj2r0_acc(xz_proj,theta_degree*pi/180,SOD,ODD,delta_T);
%add new particles
r_expression1=@(t) [4*sin(t+0.1),4*cos(t-0.3),10*t];
[xz_proj1,real_positions1] = generateTestPositions(r_expression1, conditions);
plot(xz_proj(:,1),xz_proj(:,2)); title('projection plane'); hold on
plot(xz_proj1(:,1),xz_proj1(:,2)); 
figure
plot3(real_positions(:,1),real_positions(:,2),real_positions(:,3))

