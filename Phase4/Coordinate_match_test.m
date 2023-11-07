NOS=30;
SOD = 15; ODD = 400; noise =0; delta_T = 1/68; %detector operates at 68FPStheta_degree=120*delta_T; %The tube turns 120degrees/sec
theta_degree=1.8; %The tube turns 120degrees/sec
theta=theta_degree*pi/180;
conditions=[noise, delta_T , NOS, theta, SOD , ODD];
r_expression=@(t)[1,2,3];
proj_CH=generateTestPositions(r_expression, conditions);

%% Read the csv files
folderPath='C:\Users\19606\OneDrive\Documents\GitHub\Fluid-tracking\Phase4\Variable_Speeds_Python\Data_3particles\';
FileName='Shot';
proj_Alaa=[];
for i=0:NOS-1
    fileName_i=[FileName,num2str(i)];
    fullPath = [folderPath, fileName_i,'.csv'];
    table_i=readtable(fullPath);
    proj_Alaa_i=(table_i{1:2,'Var1'})';
    proj_Alaa_i=offset(proj_Alaa_i);
    proj_Alaa=[proj_Alaa;proj_Alaa_i];
end
plot(proj_CH(:,1),proj_CH(:,2));hold on; plot(proj_Alaa(:,1),proj_Alaa(:,2))
legend('CyrusHubert','Alaa');xlabel('x_i');ylabel('y_i');


function [corrected_proj]=offset(proj)
    corrected_proj=proj;
    corrected_proj(1)=proj(1)-0;

end


