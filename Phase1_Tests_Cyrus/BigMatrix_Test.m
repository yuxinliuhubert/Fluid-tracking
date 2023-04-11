format default
xz_proj=[0.164444774551499   0.504801906568651;
   0.336777729201750   0.536135423005341;
   0.432301826101748   0.587253944887133];%theta =pi/6
theta =pi/6;
bb=proj2r0(xz_proj,theta,1,1)
function [b]= proj2r0(proj,theta,SRD,RDD) %takes in any N (even #) by 2 of projection coordinates and the angle 
   %theta btw the shots and outputs the predicted [x0, y0, z0] sets in coordinates when theta=0 (original
   %position before camera machine starts rotating)
   %r0=[x0_1, y0_1,z0_1; x0_2,y0_2,z_02;.....]
   NOS=height(proj); SDD=(SRD+RDD);
   A=zeros(4*NOS-2,1+2*NOS); A(:)=2.0000; b=zeros(height(A),1); zi_1=proj(1,2); xi_1=proj(1,1);
   A(1:2,1:3)=[1,0,-zi_1/SDD;0, -1, xi_1/SDD];
   b(1:2)=[zi_1*SRD/SDD,-xi_1*SRD/SDD];
       %b=[0;0;-xi_2*SRD/(SRD+RDD);-xi_1*SRD/(SRD+RDD);zi_1*SRD/(SRD+RDD);zi_2*SRD/(SRD+RDD)];
       for j = 2:(NOS) %j represents (number of rotations-1),  so j=1 means not rotated and j=2 means one rotation..
           %for each increased number of shots there are 2 new variables introduced
           xi_j=proj(j,1);   zi_j=proj(j,2); 
           A(4*j-5,1)=1; A( 4*j-5,2*j+1)=-zi_j/SDD; b(4*j-5)=zi_j*SRD/SDD;%for z0 magnification eq
           A(4*j-4,2*j)=-1;A(4*j-4,2*j+1)=xi_j/SDD; b(4*j-4)=-xi_j*SRD/SDD;%for x_0 magnification eq
           %the following 2 rows of codes are for equations transformation extracted from transformation
           A(4*j-3,2*j-2:2*j)=[cos(theta),sin(theta),-1]; 
          A(4*j-2,2*j-2:2*j+1)=[-sin(theta),cos(theta),0,-1]; 
          
            
            %x=(A\b);
            
          
       end
   function [r2]=T_2d(r1,alpha)
   r2=[cos(-alpha) -sin(-alpha); sin(-alpha) cos(-alpha)]*r1;
   end
end
   
