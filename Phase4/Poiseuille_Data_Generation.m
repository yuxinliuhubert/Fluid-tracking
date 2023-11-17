%% Setup parameters
format long
SOD = 15; ODD = 400;NOS=180; noise =0.01; delta_T = 1/68; %detector operates at 68FPS
theta_degree=10; 
theta=theta_degree*pi/180;
conditions=[noise, delta_T , NOS, theta, SOD , ODD];
%Fluids parameters
R=6.35/2; mu=1;dpdx=-1; %The radius of the test tube, viscocitym and pressure gradient
%% Define starting xy positions and generate projections
xy=[0,0; -0.5,0.5;-1,1 ; 1.5,1.5 ; 2,-2]; %The initial xy positions of different particles concatenated vertically
NOP=size(xy,1);%NOP stands for number of particles
xz_proj=zeros(NOS,2*NOP);real_positions=zeros(NOS,3*NOP);
for i=1:NOP
    x=xy(i,1);y=xy(i,2);
    r=sqrt(x^2+y^2);
    w=(-1/(4*mu))*dpdx*(R^2-r^2);
    r_expression=@(t) [x,y,w*t];
    [xz_proj(:,(2*i-1):2*i), real_positions(:,(3*i-2):3*i) ]=generateTestPositions(r_expression, conditions);

end
% save the data with right sequence as CSV files
folder = 'C:\Users\Crysis\Documents\Github\Fluid-tracking\Phase4\Poiseuille_Flow_Data';
for i=1:NOS
    filename_proj=[folder,'\Sorted\NOS',num2str(i-1),'.csv'];
    data_proj=xz_proj(i,:)';
    filename_pos=[folder,'\Real_positions\Pos_NOS',num2str(i-1),'.csv'];
    data_pos=real_positions(i,:)';
    csvwrite(filename_proj, data_proj);
    csvwrite(filename_pos, data_pos);
end

%% Scramble data to make unsorted
%Initializing the final scrambled matrix that has all shots and particles
scrabled_proj=zeros(NOS,NOP*2); 
for i=1:NOS
shuffle=randperm(NOP);
%this is scrambling projections for a single shot
scrambled_proj_i=zeros(1,NOP*2);
    for j=1:NOP
        scrambled_proj_i(j*2-1:j*2)=xz_proj(i,(shuffle(j)*2-1):shuffle(j)*2);
    end
    scrabled_proj(i,:)=scrambled_proj_i;
    data_proj=scrambled_proj_i';
    filename_proj=[folder,'\Unsorted\NOS',num2str(i-1),'.csv'];
    csvwrite(filename_proj, data_proj);
end