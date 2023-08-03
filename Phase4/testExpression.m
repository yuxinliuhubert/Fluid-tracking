function v = testExpression(index,T)
h = 5; % height
r = 1; % radius
 omega = 2*pi/(T/3);
 A = 1;  % initial amplitude of 1 m
switch index

    case 0
        v=@(t)[0.09*sin(t), 0.09*cos(t),0.1];


    case 1
%     r = 0.5; % radius
    h = 5; % height
    omega = 2*pi/(T/3); % angular velocity for circular motion

    v = @(t) ((t <= T/3).*([r/(T/3); 0; 0]) + ...  % First third - linear in x
     (t > T/3 & t <= 2*T/3).*([-r*omega*sin(omega*(t-T/3)); r*omega*cos(omega*(t-T/3)); 0]) + ... % Middle third - circular path in xy-plane
     (t > 2*T/3).*([-r*omega*sin(omega*(t-2*T/3)); r*omega*cos(omega*(t-2*T/3)); h/(T/3)]))'; % Final third - circular path in xy-plane and linear in z

         
    case 2
%         r = 0.5; % radius
    h = 5; % height
    omega = 2*pi/(T/3); % angular velocity for circular motion

        v = @(t) ((t <= T/3).*([-r*omega*sin(omega*t); r*omega*cos(omega*t); 0]) + ...  % First third - circular path in xy-plane
     (t > T/3 & t <= 2*T/3).*([0; 0; h/(T/3)]) + ... % Middle third - linear in z
     (t > 2*T/3).*([-r*omega*sin(omega*(t-2*T/3)); r*omega*cos(omega*(t-2*T/3)); 0]))'; % Final third - circular path in xy-plane





    case 3
            a = -0.4/(T^2/36);
             b = 0.4/T;
     
    omega = 2*pi/(T/3); % angular velocity for circular motion


        v = @(t) ((t <= T/3).*([a*t.^2; b*t; 0]) + ...  % First third - parabolic in x and linear in y
     (t > T/3 & t <= 2*T/3).*([-r*omega*sin(omega*(t-T/3)); r*omega*cos(omega*(t-T/3)); 0]) + ... % Middle third - circular path in xy-plane
     (t > 2*T/3).*([0; 0; h/(T/3)]))'; % Final third - linear in z

    case 4
            A = 1;  % initial amplitude of 1 m
%     b = 0.1;  % chosen for reasonable damping within T/3
    omega = 2*pi/(T/3);
%     r = 0.5; % radius
    v = @(t) ((t <= T/3).*([A*sin(omega*t); 0; 0]) + ...  % First third - sinusoidal in x
     (t > T/3 & t <= 2*T/3).*([0; r/(T/3); 0]) + ... % Middle third - linear in y
     (t > 2*T/3).*([-r*omega*sin(omega*(t-2*T/3)); r*omega*cos(omega*(t-2*T/3)); h/(T/3)]))'; % Final third - circular path in xy-plane and linear in z

    case 5 % damped harmonic motion
        b = 0.1;  % chosen for reasonable damping within T/3
        v = @(t) ((t <= T/3).*([A*exp(-b*t).*sin(omega*t); 0; 0]) + ...  % First third - damped sinusoidal in x
     (t > T/3 & t <= 2*T/3).*([0; r/(T/3); 0]) + ... % Middle third - linear in y
     (t > 2*T/3).*([-r*omega*sin(omega*(t-2*T/3)); r*omega*cos(omega*(t-2*T/3)); h/(T/3)]))'; % Final third - circular path in xy-plane and linear in z






end
end