function [percent_error]=Consecutive(true_position_3d,noise,NOS, theta_degree)
%% Data Input
format long
SRD = 1; % m, Source-Reference Distance
RDD = 1; % m, Reference-Detector (screen) Distance
D = 1e-4; % m, diameter of the fluid element, assuming 3D sphere/2D circle
delta_T = 1; % s, time between shots
a = [0 0]; % m/s^2, fluid element acceleration
% % True location of particle; could be negative number
x0 = true_position_3d(1); y0 = true_position_3d(2); z0=true_position_3d(3);  % in m, relative to reference

%% Calculation part 1: projection generation
theta = theta_degree/180*pi;
% % First shot: without rotation
r_now=[x0; y0;z0];
M_p = (SRD+RDD)/(SRD+r_now(2));
x_proj=[M_p*r_now(1)+randn*noise/2];
z_proj=[M_p*r_now(3)+randn*noise/2];
% % Second shot: with theta degrees rotation each time
for i = 1:NOS-1
    r_now=T(r_now,theta); %true location in the newly rotated axes
    M_p = (SRD+RDD)/(SRD+r_now(2)); % magnification of particle
    x_proj=[x_proj;M_p*r_now(1)+randn*noise/2];
    z_proj=[z_proj;M_p*r_now(3)+randn*noise/2];
end
%x_proj and z_proj are NOS by 1 vectors with all the x-coordinates projection on the screen
xz_proj=[x_proj, z_proj];
%% Calculation part2: getting the predicted values 
xyz_predicted=proj2r0(xz_proj,theta,SRD,RDD); 
avg_predicted=mean(xyz_predicted);
percent_error=abs((avg_predicted(1)-x0)/x0*100);

 %% function for rotating coordinate axes by a given angle alpha 
   function [r2]=T(r1,alpha) %r1=[x1;y1]and r2=[x2;y2] describe the same physical pt but in 2 coordinate systems with
   %theta degress rotation CCW
   r2=[cos(-alpha) -sin(-alpha) 0; sin(-alpha) cos(-alpha) 0;0 0 1]*r1;
   %note that we are rotating the coordinate axes, not the location vector of the particle
   end
  
   %% 
   function [r0]= proj2r0(proj,theta,SRD,RDD) %takes in any N (even #) by 2 of projection coordinates and the angle 
   %theta btw the shots and outputs the predicted [x0, y0, z0] sets in coordinates when theta=0 (original
   %position before camera machine starts rotating)
   %r0=[x0_1, y0_1,z0_1; x0_2,y0_2,z_02;.....]

   r0=zeros(height(proj)-1,3); %initialize the predicted positions
       for i = 1:(height(proj)-1)
           xi_1=proj(i,1); xi_2=proj(i+1,1); alpha=i*theta;
           A=[cos(theta), sin(theta), -1, 0; 
                 -sin(theta), cos(theta), 0, -1; 
                 -1, xi_1/(SRD+RDD), 0, 0,;
                 0, 0, -1, xi_2/(SRD+RDD)];
            b=[0;0;-xi_1*SRD/(SRD+RDD);-xi_2*SRD/(SRD+RDD);];
            x=(A\b);
            r0(i,1:2)=(T_2d(x(3:4), -alpha))';
            r0(i,3)=proj(i+1,2)*(SRD+x(4))/(RDD+SRD);
          
       end
   function [r2]=T_2d(r1,alpha)
   r2=[cos(-alpha) -sin(-alpha); sin(-alpha) cos(-alpha)]*r1;
   end
   end
   

end
