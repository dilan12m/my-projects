clear;
close all;
clc;
%% root locus tecchnique, autonmation
%   USER DEFINED AREA

fs  = 5; % sampling frequency in Hz
T = 1 /fs; %sampling period in seconds
G = tf(0.01758*[1, 0.876],poly([1, 0.6703]),T);% plant tf
zstar = 0.5158 + 0.4281i; % desired closed loop pole

%%%%%%%%%

% using the angle criterion
angleG = angle(evalfr(G,zstar)); %angle of plant evaluated at zstar
angleD = pi - angleG; % necessary angle of controller, D

% making range of possible andles of zero
numControllers =10;
%minAnglez1 = angleD;
%maxAnglez1 = pi;

minAnglez1 = 1.2;
maxAnglez1 = 1.9;


anglez1 = linspace(minAnglez1+10*pi/180, maxAnglez1-10*pi/180, numControllers);
anglep1 = anglez1 - angleD;

% let's determine the controller
for n = 1:numControllers
    %calculating locating of zero
    z1(n) = real(zstar) - imag(zstar)/ tan(anglez1(n));
    
    %calculating locating of pole
    p1(n) = real(zstar) - imag(zstar)/ tan(anglep1(n));
    
    % calculating controller w/o the gain parameter, K
    %%%%%%%%%%%   USER CHANGE  %%%%%%%%%%%
    D_no_k(n) = tf([1, -z1(n)],[1, -p1(n)],T);
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
xlabel('angle(z_1)'), ylabel('real part of colsed loop poles'), grid minor

% keyboard

% step response (helpful for comparing controllers)
ts = 2;
t = (0:T:4*ts)';
y = step(G_CL(5),t);
u = step(G_u_CL(5),t);

figure 
subplot(211), plot(t,u,'k.', 'markersize', 16), grid minor
xlabel('t [s]'), ylabel('control signal')

subplot(212), plot(t,y,'k.',t, ones(size(t)),'r-', 'markersize', 16), 
grid minor, xlabel('t [s]'), ylabel('output')

%ramp response
y_ramp = step(G_CL_ramp(5),t);
u_ramp = step(G_u_CL_ramp(5),t);

figure
subplot(211), plot(t,u_ramp,'k.', 'markersize', 16), grid minor
xlabel('t [s]'), ylabel('control signal')
subplot(212)
plot(t,y_ramp,'k.',t,t,'r-', 'markersize', 16), 
grid minor, xlabel('t [s]'), ylabel('output')











