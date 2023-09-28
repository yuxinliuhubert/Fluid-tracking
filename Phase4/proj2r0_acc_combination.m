function [r0]= proj2r0_acc_combination(proj,theta,SRD,RDD,delta_T) %takes in any N (even #) by 2 of projection coordinates and the angle 
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
       
       IoR=2*NOS+1; %IoR is for tracking the index of unfilled rows of the big matrix A
      for k = 2:(NOS)
          A( IoR:IoR+2*(k-1)-1, 2*k : 2*k+1 )=repmat([-1 0; 0 -1],k-1,1);
          
          for l=1:k-1
              A(IoR+2*(l-1) : IoR+2*(l-1)+1 , 2*(k-1)-2*(l-1):2*(k-1)-2*(l-1)+1)=[cos(theta*l) sin(theta*l); -sin(theta*l) cos(theta*l)];
          end
          IoR=IoR+2*(k-1);
      end
            
       
%% Now, we expand the number of columns to incorporate new variables: u, v , w, a_x, a_y, a_z
A=[A,zeros(height(A),6)];
new_col_num=size(A,2);%new_col_num is the number of columns after adding the new variables 
%u, v , w, a_x, a_y, a_z are at the last 6 columns

%% The following is for coefficients related to V and a to magnification equations 
     for j = 1:(NOS)  
           %w is at new_col_num (the last column)
           A( 2*j-1,new_col_num-3)=delta_T*(j-1);
           A( 2*j-1,new_col_num)=0.5*(delta_T*(j-1))^2;%for w term 
           
     end
     % The following is for adding equations related to u, v, a_x, a_y to
     % transformation equations
     u_ind=new_col_num-5; v_ind=u_ind+1;w_ind=u_ind+2;ax_ind=u_ind+3;ay_ind=u_ind+4;az_ind=u_ind+5;
     IoR=2*NOS+1; %IoR is for tracking the index of unfilled rows of the big matrix A
     
     for k=2:NOS  
         
         NoI=k-1; %NoT stands for the number of intervals, each with delta_T and theta, between the shot being
         % transformed and the first shot
              A(IoR, new_col_num-5)=cos(theta*NoI)*delta_T*(NoI);
              A(IoR+1, new_col_num-5)=-sin(theta*NoI)*delta_T*(NoI);%for u
              A(IoR, new_col_num-4)=sin(theta*NoI)*delta_T*(NoI);
              A(IoR+1, new_col_num-4)=cos(theta*NoI)*delta_T*(NoI);%for v
              A(IoR, new_col_num-2)=0.5*cos(theta*NoI)*(delta_T*NoI)^2;
              A(IoR+1, new_col_num-2)=-sin(theta*NoI)*0.5*(delta_T*NoI)^2;%for ax
              A(IoR, new_col_num-1)=0.5*sin(theta*NoI)*(delta_T*NoI)^2;
              A(IoR+1, new_col_num-1)=0.5*cos(theta*NoI)*(delta_T*NoI)^2;%for ay

             
            
         
          IoR=IoR+2;
     end
     IoR=2*NOS+1;
     for k = 2:(NOS)
          for l=1:k-1
            NoI=k-l;
              A(IoR, x_ind)=cos(theta*NoI)*delta_T*(NoI);
              A(IoR+1, new_col_num-5)=-sin(theta*NoI)*delta_T*(NoI);%for u
              A(IoR, new_col_num-4)=sin(theta*NoI)*delta_T*(NoI);
              A(IoR+1, new_col_num-4)=cos(theta*NoI)*delta_T*(NoI);%for v
              A(IoR, new_col_num-2)=0.5*cos(theta*NoI)*(delta_T*NoI)^2;
              A(IoR+1, new_col_num-2)=-sin(theta*NoI)*0.5*(delta_T*NoI)^2;%for ax
              A(IoR, new_col_num-1)=0.5*sin(theta*NoI)*(delta_T*NoI)^2;
              A(IoR+1, new_col_num-1)=0.5*cos(theta*NoI)*(delta_T*NoI)^2;%for ay

          end
          IoR=IoR+2*(k-1);
      end

     x=(A\b);
r0=[x(2),x(3),x(1),x(new_col_num-5),x(new_col_num-4),x(new_col_num-3),x(new_col_num-2),x(new_col_num-1),x(new_col_num)];
   end