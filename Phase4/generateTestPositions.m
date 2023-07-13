function [xz_proj,real_positions] = generateTestPositions(vel_expression,initial_position_3d,conditions)

%% Calculation part 1: projection generation

[noise,delta_T,NOS,theta_degree,SRD,RDD] = deal(conditions(1),conditions(2), ...
    conditions(3),conditions(4),conditions(6),conditions(7));

theta = theta_degree/180*pi;
% % First shot: without rotation
r0_0=initial_position_3d;
real_positions=zeros(NOS,3);
real_positions(1,:)= initial_position_3d;
M_p = (SRD+RDD)/(SRD+r0_0(2));
x_proj=zeros(NOS,1); z_proj=zeros(NOS,1);
x_proj(1)=M_p*r0_0(1)+randn*noise/2;
z_proj(1)=M_p*r0_0(3)+randn*noise/2;
% % Second shot: with theta degrees rotation each time
for k = 1:NOS-1
    r0_k=r0_0+integral(vel_expression,0,delta_T*k,'ArrayValued', true); % true location in the original frame of reference
    real_positions(k+1,:) =  r0_k;
    r_now=T(r0_k',theta*k)';
    M_p = (SRD+RDD)/(SRD+r_now(2)); % magnification of particle
    x_proj(k+1)=M_p*r_now(1)+randn*noise/2;
    z_proj(k+1)=M_p*r_now(3)+randn*noise/2;
end
% x_proj and z_proj are NOS by 1 vectors with all the x-coordinates projection on the screen
xz_proj=[x_proj, z_proj];

end
