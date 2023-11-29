NOS=146;
SOD = 38; ODD = 500-38; noise =0; delta_T = 1/68; %detector operates at 68FPStheta_degree=120*delta_T; %The tube turns 120degrees/sec
theta_degree=1.8; %The tube turns 120degrees/sec
theta=theta_degree*pi/180;
conditions=[noise, delta_T , NOS, theta, SOD , ODD];
r_expression=@(t)[1,2,3];
r_expression1=@(t)[0,0,0]; 
proj_CH=generateTestPositions(r_expression, conditions);

%% Read the csv files
folderPath='C:\Users\19606\OneDrive\Documents\GitHub\Fluid-tracking\Phase7_Python\5_Nov_23\';
FileName='Shot';
proj_Alaa=[];
for i=0:NOS-1
    fileName_i=[FileName,num2str(i)];
    fullPath = [folderPath, fileName_i,'.csv'];
    table_i=readtable(fullPath);
    proj_Alaa_i=(table_i{1:2,'Var1'})';
proj_Alaa_i=offset(proj_Alaa_i,i-1);
    proj_Alaa=[proj_Alaa;proj_Alaa_i];
end
colors=1:1:NOS;colors=colors';
scatter(proj_CH(:,1),proj_CH(:,2),30,colors);hold on; scatter(proj_Alaa(:,1),proj_Alaa(:,2),30,colors,'filled')
legend('CyrusHubert','Alaa');xlabel('x_i');ylabel('y_i');
% Add a colorbar to show the color scale
colorbar;

% Apply a colormap
colormap(jet);

function [corrected_proj]=offset(proj,theta_count)
    corrected_proj=proj;
    corrected_proj(1)=proj(1)*0.172-243.5*0.172;
    %(proj(1)-243.5*0.175)*0.175
    corrected_proj(2)=proj(2)*0.172-97.5*0.172;

end

