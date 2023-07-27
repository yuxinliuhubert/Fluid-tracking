function v = testExpression(index,T)
switch index

    case 0
        v=@(t)[0.9*sin(t), 0.9*cos(t),1];


    case 1
        a = 9;

        v = @(t) ((t <= T/3).*([a*t; a*t; a*t] ) + ...  % First third - linear
         (t > T/3 & t <= 2*T/3).*([a*cos(2*pi*(t-T/3)/(T/3)); a*sin(2*pi*(t-T/3)/(T/3)); a*(t-T/3)/(T/3)]) + ... % Middle third - circular path in xy-plane
         (t > 2*T/3).*([a*cos(2*pi); a*sin(2*pi); a*(t-2*T/3)/(T/3)] + [a*cos(2*pi*(t-2*T/3)/(T/3)); a*sin(2*pi*(t-2*T/3)/(T/3)); a]))'; % Final third - circular path in xy-plane starting from the end of the previous path
         
    case 2
        a = 5;
        v = @(t) (t <= T/3).*(a*t.^2) + ...  % First third - quadratic
         (t > T/3 & t <= 2*T/3).*(a*sin(2*pi*(t-T/3)/(T/3))) + ...  % Middle third - sine wave
         (t > 2*T/3).*(a*log(t-2*T/3+1));  % Final third - logarithmic


%     case 3
%         v = @(t) (t<delta_T*round(NOS/2)).* [2;1*t+3*t^3;2] + ...
%         (t>=delta_T*round(NOS/2) & t < delta_T*round(NOS*0.8)) .* [3*t^3+7;2;1*t^2+2.2] + ...
%         (t>=delta_T*round(NOS*0.8)).* [2-4*t^2;1*t;2];




end
end