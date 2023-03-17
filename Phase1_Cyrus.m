% M_particle = (SRD+RDD)/(SRD+y0); % magnification of particle
%% Data Input
SRD = 1; % m, Source-Reference Distance
RDD = 1; % m, Reference-Detector (screen) Distance
theta_degree = 30; % clock-wise degree, what is the angle of camera rotation before each shot
NOS = 10; % number of shots
D = 1e-4; % m, diameter of the fluid element, assuming 3D sphere/2D circle
delta_T = 1; % s, time between shots
a = [0 0]; % m/s^2, fluid element acceleration
% % True location of particle; could be negative number
x0 = 0.1; % m, relative to reference
y0 = 0.1; % m
noise=5e-4; %the maximum possible deviation of the projection measurement from the perfect measurement 

%% Calculation part 1: projection generation
theta = theta_degree/180*pi;
%rotated basis  to original basis. For every shot taken, A multiplies by once. 
% projected x input; could be used to create kalman filter with rotation
% matrix to minimize the effect of noise. This version assumes no noise
% from measurement
% M >1
x_proj = [ % [particle left projection, particle right side projection]
    ];
r_now=[x0; y0];
M_p = (SRD+RDD)/(SRD+r_now(2));
x_proj=[M_p*r_now(1)+randn*noise];
for i = 1:NOS-1
  
    r_now=T(r_now,theta); %true location in the newly rotated axes
    M_p = (SRD+RDD)/(SRD+r_now(2)); % magnification of particle
    x_proj=[x_proj;M_p*r_now(1)+randn*noise];%x_proj=[x_proj; (r_now(1)-D/2)*M_p+randn*noise/2, (r_now(2)+D/2)*M_p+randn*noise/2];
end
%x_proj_L=x_proj(:,1) ; x_proj_R=x_proj(:,2) ; 
%%Calculation part 2: location prediction based on projection with noise
% r0_L_predicted=proj2r0(x_proj_L,theta);
% r0_R_predicted=proj2r0(x_proj_R,theta);
%r0_predicted=(r0_L_predicted+r0_R_predicted)/2
r0_predicted=proj2r0(x_proj,theta)

   %% function for transforming rotating coordinate axes of a given angle alpha and r1=[x1;y1]
   function [r2]=T(r1,alpha)
   r2=[cos(-alpha) -sin(-alpha); sin(-alpha) cos(-alpha)]*r1;
   end
   
   function [r0]= proj2r0(proj,theta) %takes in any array of projection coordinates and the angle theta btw the shots and 
   %outputs the predicted x0 y0 sets in coordinates when theta=0 (original
   %position) r0=[x0_1, y0_1; x0_2,y0_2;.....]
   SRD = 1; 
RDD = 1;
   r0=zeros(length(proj)/2,2);
       for i = 1:(length(proj)/2)
           xi_1=proj(2*i-1); xi_2=proj(2*i); alpha=2*i*theta;
           A=[cos(theta), sin(theta), -1, 0; 
                 -sin(theta), cos(theta), 0, -1; 
                 -1, xi_1/(SRD+RDD), 0, 0,;
                 0, 0, -1, xi_2/(SRD+RDD)];
            b=[0;0;-xi_1*SRD/(SRD+RDD);-xi_2*SRD/(SRD+RDD);];
            x=(A\b);
            x=(T(x(3:4), -alpha+theta))';
            r0(i,:)=x;
          
       end
   function [r2]=T(r1,alpha)
   r2=[cos(-alpha) -sin(-alpha); sin(-alpha) cos(-alpha)]*r1;
   end
   end
