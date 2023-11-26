function [xz_proj,real_positions] = generateTestPositions(r_expression, conditions)



[noise, delta_T , NOS, theta, SOD , ODD,hold] = deal(conditions(1),conditions(2), conditions(3),conditions(4),conditions(5),conditions(6),conditions(7));


%% Initialize the variables

real_positions=zeros(NOS,3);
x_proj=zeros(NOS,1); z_proj=zeros(NOS,1);

%% Loop through all the shots
    hold=abs(round(hold));
    for k = 1:NOS
        if k> hold
            r0_k=r_expression((k-1)*delta_T); % true location in the original frame of reference
            real_positions(k,:) =  r0_k;
            r_rotated=T(r0_k',theta*(k-hold))';
            M_p = (SOD+ODD)/(SOD+r_rotated(2)); % magnification of particle
            x_proj(k)=M_p*r_rotated(1)+randn*noise/2;
            z_proj(k)=M_p*r_rotated(3)+randn*noise/2;
        else
            r0_k=r_expression((k-1)*delta_T); % true location in the original frame of reference
            real_positions(k,:) =  r0_k;
            M_p = (SOD+ODD)/(SOD+r0_k(2)); % magnification of particle
            x_proj(k)=M_p*r0_k(1)+randn*noise/2;
            z_proj(k)=M_p*r0_k(3)+randn*noise/2;
        end
    end
    % x_proj and z_proj are NOS by 1 vectors with all the x-coordinates projection on the screen
    xz_proj=[x_proj, z_proj];
    
end
