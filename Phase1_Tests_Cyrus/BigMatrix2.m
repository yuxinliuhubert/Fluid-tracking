function [percent_error]=BigMatrix2(true_position_3d,noise,NOS, theta_degree)
%% Data Input
format long
SRD = 1; % m, Source-Reference Distance
RDD = 1; % m, Reference-Detector (screen) Distance
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
percent_error=abs((xyz_predicted(1)-x0)/x0*100);


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
   NOS=height(proj); SDD=(SRD+RDD);
   row_number_A=round(2*NOS+ 2*(factorial(NOS)/(factorial(NOS-2)*2))  ,  0  ); 
   col_number_A=round(1+2*NOS, 0); %round function is for avoiding precision error
   A=zeros(   row_number_A,  col_number_A); 
   b=zeros(height(A),1);
       for j = 1:(NOS)  %This for loop is for constructing the equations arising from magnification alone
           %for each increased number of shots there are 2 new variables
           %introduced and 2 equations
           xi_j=proj(j,1);   zi_j=proj(j,2); 
           A(2*j-1,1)=1; A( 2*j-1,2*j+1)=-zi_j/SDD; b(2*j-1)=zi_j*SRD/SDD;%for z0 magnification eq
           A(2*j,2*j)=-1;A(2*j,2*j+1)=xi_j/SDD; b(2*j)=-xi_j*SRD/SDD;%for x_0 magnification eq
           %the following 2 rows of codes are for equations transformation extracted from transformation
       end
       
       x=2*NOS+1; %x is for tracking the index of unfilled rows of the big matrix A
      for k = 2:(NOS)
          A( x:x+2*(k-1)-1, 2*k : 2*k+1 )=repmat([-1 0; 0 -1],k-1,1);
          
          for l=1:k-1
              A(x+2*(l-1) : x+2*(l-1)+1 , 2*(k-1)-2*(l-1):2*(k-1)-2*(l-1)+1)=[cos(theta*l) sin(theta*l); -sin(theta*l) cos(theta*l)];
          end
          x=x+2*(k-1);
      end
            x=(A\b);
            r0=[x(2),x(3),x(1)];
          
       
 
   end
   

end
