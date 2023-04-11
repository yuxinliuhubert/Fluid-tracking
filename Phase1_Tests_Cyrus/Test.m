true_position_3d=[0.1,0.1,0.1];NOS=30; theta_degree=30;
Runs=1000; %number of runs for each value of NOS
noise_ratio=1; %noise_ratio is the ratio of position value to the error for each projection 
% for 68% of the time (normal distribution); note that this ratio is not in percent so 0.1 means 10%
noise=mean(true_position_3d)*noise_ratio;


%% Test for the original method: 1,2  3,4  5,6.....
%NOS for the original method has to be a even number
Original_Result_PercentError=zeros(NOS/2,1);

for i =2:2:NOS
     Error_ThisNOS=zeros(Runs,1);
    for j=1:Runs
    Error_ThisNOS(j)=Original(true_position_3d,noise,NOS, theta_degree);
    end
    Original_Result_PercentError(i/2)=mean(Error_ThisNOS);
end

%% Test for consecutive number average method: 1,2  2,3 3,4...
Consecutive_Result_PercentError=zeros(NOS-1,1);
for i =2:NOS
    Error_ThisNOS=zeros(Runs,1);
    for j=1:Runs
    Error_ThisNOS(j)=Consecutive(true_position_3d,noise,NOS, theta_degree);
    end
    Consecutive_Result_PercentError(i-1)=mean(Error_ThisNOS);
end

%% Test for original method but the matrix is bigger for including equations introduced by the added z-axis
OriginalBigger_Result_PercentError=zeros(NOS/2,1);

for i =2:2:NOS
    Error_ThisNOS=zeros(Runs,1);
    for j=1:Runs
    Error_ThisNOS(j)=Original_BiggerMatrix(true_position_3d,noise,NOS, theta_degree);
    end
    OriginalBigger_Result_PercentError(i/2)=mean(Error_ThisNOS);
end

%% Test for the big matrix method
BigMatrix_Result_PercentError=zeros(NOS/2,1);

for i =2:NOS
    Error_ThisNOS=zeros(Runs,1);
    for j=1:Runs
    Error_ThisNOS(j)=BigMatrix2(true_position_3d,noise,NOS, theta_degree);
    end
    BigMatrix_Result_PercentError(i-1)=mean(Error_ThisNOS);
end


%% Plotting
figure
plot(linspace(2,NOS,NOS/2),Original_Result_PercentError)
hold on
plot(linspace(2,NOS,NOS-1),Consecutive_Result_PercentError)
%plot(linspace(2,NOS,NOS/2),OriginalBigger_Result_PercentError)
plot(linspace(2,NOS,NOS-1),BigMatrix_Result_PercentError)
title(['Rotation degree=', num2str(theta_degree),' Nominal Error=',num2str(noise)] )
xlabel('Number of Shots')
ylabel('percent error')
legend('Original','Consecutive','BigMatrix')
%legend('Original','Consecutive','OriginalBig','BigMatrix')

