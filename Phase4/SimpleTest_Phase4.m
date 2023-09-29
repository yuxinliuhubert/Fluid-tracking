SOD = 1; RDD = 1;noise = 0; delta_T = 0.1; NOS=5;theta_degree=10;
theta=theta_degree*pi/180;
conditions=[noise,delta_T,NOS,theta_degree,SOD,RDD];
r_expression=@(t) [0.4+t,0,0.05*t^2+0.3*t+1];
[xz_proj,real_positions] = generateTestPositions(r_expression, conditions);
correct_ans=proj2r0_acc(xz_proj,theta_degree*pi/180,SOD,RDD,delta_T);
