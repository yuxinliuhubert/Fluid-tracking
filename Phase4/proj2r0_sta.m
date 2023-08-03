   function [r0]= proj2r0_sta(proj,theta,SRD,RDD) %takes in any N (even #) by 2 of projection coordinates and the angle 
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