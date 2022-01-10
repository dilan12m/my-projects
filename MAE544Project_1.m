clear;
close all;
clc;
%% root locus tecchnique, autonmation
%   USER DEFINED AREA

fs = 20; % sampling frequency in Hz
T  = 1 /fs; %sampling period in seconds
G  = tf(1.796,[1, -0.5134],T); % plant tf
kp = tf(1,[1, -1],T);

% desired closed loop pole
%zstar = .62492202 + .37805122i; 
%zstar = .7325540 + .186136i;
%zstar = .648 + .311i;
%zstar = .731475021 + .1858684025i;
zstar = .65855 + .28162i;


%zstar = .661 + .277i;


%%%%%%%%%

% using the angle criterion
angleG = angle(evalfr(G,zstar)); %angle of plant evaluated at zstar
angleD = pi - angleG; % necessary angle of controller, D
p2 = angle(evalfr(kp,zstar));

% making range of possible angles of zero
numControllers = 40;
minAnglez1 = angleD;
maxAnglez1 = pi;

anglez1 = linspace(minAnglez1+10*pi/180, maxAnglez1-10*pi/180, numControllers);
anglep1 = anglez1 - angleD + p2;

% let's determine the controller
for n = 1:numControllers
    %calculating location of zero
    z1(n) = real(zstar) - imag(zstar)/ tan(anglez1(n));
    
    %calculating location of pole
    p1(n) = real(zstar) - imag(zstar)/ tan(anglep1(n));
    
    
    % calculating controller w/o the gain parameter, K
    %%%%%%%%%%%   USER CHANGE  %%%%%%%%%%%
    %D_no_k(n) = tf([1, -z1(n)],[1,(-1-p1(n)),p1(n)],T);  
    D_no_k(n) = tf([1, -z1(n)],[1, -p1(n)],T)*tf(1,[1, -1],T);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %using magnitude carterion to get gain, K
    K(n) = 1/abs(evalfr(G*D_no_k(n), zstar));
    
    %put together the controller
    D(n) = K(n)*D_no_k(n);
    
    %calculate the closed-loop poles
    G_CL(n) = feedback(G*D(n),1);
    CL_poles(:,n) = pole (G_CL(n));
    
    %closed-loop control signal transfer function
    G_u_CL(n) = D(n) / (1+G*D(n));
    G_u_CL_ramp(n) = G_u_CL(n)* tf(T,[1, -1],T);
    
    %closed-loop transfer function for calculating ramp response
    G_CL_ramp(n) = G_CL(n)* tf(T,[1, -1],T);
    
    % Debug
    angle(evalfr(G*D(n),zstar));
    abs((evalfr(G*D(n),zstar)));
    
end

% visualizing closed-loop poles
CL_poles_real = real( CL_poles);
CL_poles_imag = imag( CL_poles);

figure 
plot (CL_poles_real,CL_poles_imag, 'rx', 'markersize', 16, 'linewidth',3)
grid minor;
zgrid, axis equal

figure 
plot(anglez1, CL_poles_real','rx', 'markersize', 16, 'linewidth',3);
xlabel('(z_1)'), ylabel('real part of colsed loop poles'), grid minor

% keyboard

i=15;
% step response (helpful for comparing controllers)
ts = 3;
t = (0:T:ts)';
y = 2*pi*step(G_CL(i),t);
u = 2*pi*step(G_u_CL(i),t);

figure 
subplot(211), plot(t,u,'k.', 'markersize', 16), grid minor
xlabel('t [s]'), ylabel('control signal')

subplot(212), plot(t,y,'k.',t, 2*pi*ones(size(t)),'r-', 'markersize', 16), 
grid minor, xlabel('t [s]'), ylabel('output')

%ramp response
y_ramp = step(G_CL_ramp(i),t);
u_ramp = step(G_u_CL_ramp(i),t);

figure
subplot(211), plot(t,u_ramp,'k.', 'markersize', 16), grid minor
xlabel('t [s]'), ylabel('control signal')
subplot(212)
plot(t,y_ramp,'k.',t,t,'r-', 'markersize', 16), 
grid minor, xlabel('t [s]'), ylabel('output')
title ('Ramp response')





















