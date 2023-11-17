function [xz_proj,real_positions] = generateTestPositions(r_expression, conditions)



[noise, delta_T , NOS, theta, SOD , ODD] = deal(conditions(1),conditions(2), conditions(3),conditions(4),conditions(5),conditions(6));


%% Initialize the variables

real_positions=zeros(NOS,3);
x_proj=zeros(NOS,1); z_proj=zeros(NOS,1);

%% Loop through all the shots
for k = 1:NOS
    r0_k=r_expression((k-1)*delta_T); % true location in the original frame of reference
    real_positions(k,:) =  r0_k;
    r_rotated=T(r0_k',theta*(k-1))';
    M_p = (SOD+ODD)/(SOD+r_rotated(2)); % magnification of particle
    x_proj(k)=M_p*r_rotated(1)+randn*noise/2;
    z_proj(k)=M_p*r_rotated(3)+randn*noise/2;
end
% x_proj and z_proj are NOS by 1 vectors with all the x-coordinates projection on the screen
xz_proj=[x_proj, z_proj];

end
