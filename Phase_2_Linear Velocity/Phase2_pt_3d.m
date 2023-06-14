%% Data Input
function results = Phase2_pt_3d(initial_position_3d, noise, delta_T, NOS,theta_degree)
format default
SRD = 1; % m, Source-Reference Distance
RDD = 1; % m, Reference-Detector (screen) Distance
% % theta_degree = 10; % clock-wise degree, what is the angle of camera rotation before each shot
% % NOS = 100; % number of shots
% % delta_T = time between shots
% % initial_position_3d = [x0;y0;z0], a column vector
% noise=5e-3; %the 68% chance deviation of the projection measurement from the perfect measurement 

method = 0; % 0 for least square, 1 for kalman

v = @(t) [0.1;0.3;0.2]; %velocity function

% Process error
delta_P_X = 1e-3; % m
delta_P_VX = 1e-3; % m/sec

% Observation error
delta_O_X= 2e-3; % m
delta_O_VX = 1.5e-4; % m
% delta_O_Y = delta_O_X;
% delta_O_VY = delta_O_VX;

%% Calculation part 1: projection generation
theta = theta_degree/180*pi;
% % First shot: without rotation
r0_0=initial_position_3d;
real_Positions = initial_position_3d;
M_p = (SRD+RDD)/(SRD+r0_0(2));
x_proj=M_p*r0_0(1)+randn*noise/2;
z_proj=M_p*r0_0(3)+randn*noise/2;
% % Second shot: with theta degrees rotation each time
for k = 1:NOS-1
    r0_k=r0_0+integral(v,0,delta_T*k,'ArrayValued', true); % true location in the original frame of reference
    real_Positions = [real_Positions; r0_0 + integral(v,0,delta_T*k,'ArrayValued', true)'];
    r_now=T(r0_k,theta*k);
    M_p = (SRD+RDD)/(SRD+r_now(2)); % magnification of particle
    x_proj=[x_proj;M_p*r_now(1)+randn*noise/2];
    z_proj=[z_proj;M_p*r_now(3)+randn*noise/2];
end
% x_proj and z_proj are NOS by 1 vectors with all the x-coordinates projection on the screen
xz_proj=[x_proj, z_proj];



%% Calculation part2: geting the measured values 


proj2r0(xz_proj,theta,SRD,RDD,delta_T)


 %% function for rotating coordinate axes by a given angle alpha 
   function [r2]=T(r1,alpha) %r1=[x1;y1;z1]and r2=[x2;y2;z2] describe the same physical pt but in 2 coordinate systems with
   % theta degress rotation CCW
   r2=[cos(-alpha) -sin(-alpha) 0; sin(-alpha) cos(-alpha) 0;0 0 1]*r1;
   % note that we are rotating the coordinate axes, not the location vector of the particle
   end
  
  
%row_number_A=row_number_A+NOS ;%which is different from phase 1. We get a new equation each shot for the change of z coordinate due to velocity

function [r0]= proj2r0(proj,theta,SRD,RDD,delta_T) %takes in any N (even #) by 2 of projection coordinates and the angle 
   %theta btw the shots and outputs the predicted [x0, y0, z0] sets in coordinates when theta=0 (original
   %position before camera machine starts rotating
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
       
%% Now, we expand the number of columns to incorporate new variables: u, v , w
A=[A,zeros(height(A),3)];
new_col_num=size(A,2); %new_col_num is the number of columns after adding the new variables 
%u, v , w are at the last 6 columns
%% The following is for coefficients related to V to magnification equations 
     for j = 1:(NOS)  
           %w is at new_col_num (the last column)
           A( 2*j-1,new_col_num)=delta_T*(j-1); %for w term 
           
     end
     % The following is for adding equations related to u, v, a_x, a_y to
     % transformation equations
     IoR=2*NOS+1; %IoR is for tracking the index of unfilled rows of the big matrix A
     for k=2:NOS  
         
          for L=1:k-1
              A(IoR+2*(L-1), new_col_num-2)=delta_T*(L-1);
              A(IoR+2*(L-1)+1, new_col_num-2)=delta_T*(L-1);%for u
              A(IoR+2*(L-1), new_col_num-1)=delta_T*(L-1);
              A(IoR+2*(L-1)+1, new_col_num-1)=delta_T*(L-1);%for v
             

          end
          IoR=IoR+2*(k-1);
     end

     x=(A\b);
     r0=[x(2),x(3),x(1),x(new_col_num-2),x(new_col_num-1),x(new_col_num)];
   end

   

end
