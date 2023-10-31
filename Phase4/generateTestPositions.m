function [xz_proj,real_positions] = generateTestPositions(r_expression, conditions)



[noise, delta_T , NOS, theta, SOD , ODD] = deal(conditions(1),conditions(2), conditions(3),conditions(4),conditions(5),conditions(6));


%% First shot: without rotation
r0_0=r_expression(0); %initial position
real_positions=zeros(NOS,3);
real_positions(1,:)= r0_0;
M_p = (SOD+ODD)/(SOD+r0_0(2));
x_proj=zeros(NOS,1); z_proj=zeros(NOS,1);
x_proj(1)=M_p*r0_0(1)+randn*noise/2;
z_proj(1)=M_p*r0_0(3)+randn*noise/2;
% % Second shot: with theta degrees rotation each time
for k = 1:NOS-1
    r0_k=r_expression(k*delta_T); % true location in the original frame of reference
    real_positions(k+1,:) =  r0_k;
    r_rotated=T(r0_k',theta*k)';
    M_p = (SOD+ODD)/(SOD+r_rotated(2)); % magnification of particle
    x_proj(k+1)=M_p*r_rotated(1)+randn*noise/2;
    z_proj(k+1)=M_p*r_rotated(3)+randn*noise/2;
end
% x_proj and z_proj are NOS by 1 vectors with all the x-coordinates projection on the screen
xz_proj=[x_proj, z_proj];

end
