function [r0]= proj2r0_acc_combination(xz_proj,theta,SOD,ODD,dt) %takes in any N (even #) by 2 of projection coordinates and the angle 
   %theta btw the shots and outputs the predicted [x0, y0, z0] sets in coordinates when theta=0 (original
   %position before camera machine starts rotating
    NOS=height(xz_proj); SDD=(SOD+ODD);
     row_number_A=round(2*NOS+ 2*(factorial(NOS)/(factorial(NOS-2)*2))  ,  0  ); 
   col_number_A=round(1+2*NOS, 0); %round function is for avoiding precision error
   A=zeros(   row_number_A,  col_number_A); 
   b=zeros(height(A),1);
%A=sym(A); %for debugging purpose
   %% This for loop is for constructing the equations arising from magnification alone
       for j = 1:(NOS)  
           %for each increased number of shots there are 2 new variables
           %introduced and 2 equations
           xi_j=xz_proj(j,1);   zi_j=xz_proj(j,2); 
           A(2*j-1,1)=1; A( 2*j-1,2*j+1)=-zi_j/SDD; b(2*j-1)=zi_j*SOD/SDD;%for z0 magnification eq
           A(2*j,2*j)=-1;A(2*j,2*j+1)=xi_j/SDD; b(2*j)=-xi_j*SOD/SDD;%for x_0 magnification eq
           %the following 2 rows of codes are for equations transformation extracted from transformation
       end
  %% This loop is for constructing equations from rotational transformation     
       IoR=2*NOS+1; %IoR is for tracking the index of unfilled rows of the big matrix A
       
      for k = 2:(NOS)
          trans_count=k-1; %trans_count is the number of transformations available to k's position
          A( IoR   :  (IoR-1)+(trans_count*2)   ,    2*k  :  2*k+1 )=repmat([-1 0; 0 -1],trans_count,1); %placing -1* identity matrix
          %for the k's position
          
          for i=1:trans_count %index i corresponds to the index for the position being transformed to k
              
              delta_theta=(theta*i); %theta_multi stands for multiples of theta, or number of theta's 
              A(IoR+2*(i-1)  :  IoR+2*(i-1)+1 ,  2*trans_count-2*(i-1)  :  2*trans_count-2*(i-1)+1)=[cos(delta_theta) sin(delta_theta); -sin(delta_theta) cos(delta_theta)];
          end
          IoR=IoR+2*(k-1);
      end
            
       
%% Now, we expand the number of columns to incorporate new variables: u, v , w, a_x, a_y, a_z
A=[A,zeros(height(A),6)];
new_col_num=size(A,2);%new_col_num is the number of columns after adding the new variables 
%u, v , w, a_x, a_y, a_z are at the last 6 columns
u_ind=new_col_num-5; v_ind=u_ind+1;w_ind=u_ind+2;ax_ind=u_ind+3;ay_ind=u_ind+4;az_ind=u_ind+5;
%% The following is for coefficients related to V and a to magnification equations 
%only z-directions needs such correction
     for j = 1:(NOS)  
           %w is at new_col_num (the last column)
           A( 2*j-1,w_ind)=dt*(j-1);%for w term
           A( 2*j-1,az_ind)=0.5*(dt*(j-1))^2; %for a_z term
           
     end
     %% The following is for adding equations related to u, v, a_x, a_y to
     % transformation equations
     
     IoR=2*NOS+1; %IoR is for tracking the index of unfilled rows of the big matrix A
     
     for k = 2:(NOS)
         trans_count=k-1;      
         T2=dt*(k-1);
         for i=1:trans_count %index i corresponds to the index for the position being transformed to k
              delta_theta=theta*i; %delta_theta is the angle between the frames being transformed;
              %it's multiples of theta, or number of theta's that's between j and k
             T1=T2-i*dt;
             theta_prime= theta*(k-i-1);%theta_prime stands for the multiple of theta's from the 1st frame to j frame. 
             % This needs to be distinguished fro theta_multi
             %to add coeffcients for u
              A(IoR, u_ind)=(cos(theta_prime)*cos(delta_theta)-sin(theta_prime)*sin(delta_theta))*(T2-T1);
              A(IoR+1, u_ind)=(-cos(theta_prime)*sin(delta_theta)-sin(theta_prime)*cos(delta_theta))*(T2-T1);
              %for v
              A(IoR, v_ind)=(sin(theta_prime)*cos(delta_theta)+cos(theta_prime)*sin(delta_theta))*(T2-T1);
              A(IoR+1, v_ind)=(-sin(theta_prime)*sin(delta_theta)+cos(theta_prime)*cos(delta_theta))*(T2-T1);
              %for ax
              A(IoR, ax_ind)=0.5*(cos(theta_prime)*cos(delta_theta)-sin(theta_prime)*sin(delta_theta))*(T2^2-T1^2);
              A(IoR+1, ax_ind)=0.5*(-cos(theta_prime)*sin(delta_theta)-sin(theta_prime)*cos(delta_theta))*(T2^2-T1^2);
              %for ay
              A(IoR, ay_ind)=0.5*(sin(theta_prime)*cos(delta_theta)+cos(theta_prime)*sin(delta_theta))*(T2^2-T1^2);
              A(IoR+1, ay_ind)=0.5*(-sin(theta_prime)*sin(delta_theta)+cos(theta_prime)*cos(delta_theta))*(T2^2-T1^2);
            IoR=IoR+2;
         end
       
     end
      



     x=(A\b);
r0=[x(2),x(3),x(1),x(u_ind),x(v_ind),x(w_ind),x(ax_ind),x(ay_ind),x(az_ind)];
   end