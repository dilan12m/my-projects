clear;
close all;
clc;
%% State and Initialize

dt = 0.1;
t  = [0:dt:100]';
m  = length(t);
% Constant Measurement Matrix
h  = [1 0 0 0; 0 0 1 0];

% Preallocating Space
xe    = zeros(m,4);
x     = zeros(m,4);
p_cov = zeros(m,4);
ym    = zeros(m,2);

% Gravational Constant
G   = 6.6742*10^-11; 
% Mass of the Earth
M   = 5.98 *10^24;
% Radius of the Earth
r_0 = 6.57*10^6;
% Angular Velocity
w   = sqrt((G*M)/r_0^3);
% Initial Conditions
x0  = [6570000;0;0;1.05*w];

%% Measurement and Process Noise

% Standard Deviation Squared (Variance) of Measurement and Process Noise
% Respectively
sigy = (0.1)^2; %% 1-> 0.1, 2-> 0.001, 3-> 0.0,  4-> 0.1
sigx = (0.1)^2; %% 1-> 0.1  2-> 0.001, 3-> 0.1,  4-> 0.0

% Measurement Noise
r  = sigy* eye(2);
% Process Noise 
q  = sigx*eye(4);

% Covariance
p0 = .1^2*eye(4);
p  = p0;
p_cov(1,:) = diag(p0)';

x (1,:) = x0';
xe(1,:) = x0';

% Main Routine
for i=1:m-1

    % Truth and Measurements
    % Using this Method instead of ode45

    f1 = dt*rungakutta(x(i,:));
    f2 = dt*rungakutta(x(i,:)+0.5*f1');
    f3 = dt*rungakutta(x(i,:)+0.5*f2');
    f4 = dt*rungakutta(x(i,:)+f3');
    x(i+1,:) = x(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');

    y  = (h*x')';
    ym = y+sqrt(sigy)*randn(1); % Measurement Noise 

    % Kalman Update
    k = p*h'*inv(h*p*h'+r);
    p = (eye(4)-k*h)*p;
    xe(i,:) = xe(i,:)+(k*(ym(i,:)'- h*xe(i,:)'))';

    % Propagation
    f1 = dt*rungakutta(xe(i,:));
    f2 = dt*rungakutta(xe(i,:)+0.5*f1');
    f3 = dt*rungakutta(xe(i,:)+0.5*f2');
    f4 = dt*rungakutta(xe(i,:)+f3');
    xe(i+1,:) = xe(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');

    % Linearized Matrix 
    fpart = [0 1 0 0;... 
             (xe(i,4)^2)+((2*(6.6742*10^-11)*(5.98 *10^24))/(xe(i,1)^3)) 0 0 2*xe(i,1)*xe(i,4);...
             0 0 0 1;...
             (2*xe(i,4)*xe(i,2))/(xe(i,1)^2) -2*xe(i,4)/xe(i,1) 0 (-2*xe(i,2))/(xe(i,1))];

    % Using c2d to get discrete time system
    phi = c2d(fpart,[6570000;0;0;1.05*w],dt);
    p   = phi*p*phi'+ q*dt;
    p_cov(i+1,:) = diag(p)';

end



% 3-Sigma Outlier
sig3 = p_cov.^(0.5)*3;

%% Other Useful metrics 

% Mean of the position error 
mean_position = mean(xe(:,1)-x(:,1))

% Mean of the attitude error
mean_attitude = mean(xe(:,3)-x(:,3))


%% Plot Results

figure (1)
plot(t,sig3(:,1),'r',t,xe(:,1)-x(:,1),'b',t,-sig3(:,1),'r')
xlabel('t'), ylabel('Position Errors')
title('Position of Earth Orbiting Satellite')

figure(2)
plot(t,sig3(:,3),'r',t,xe(:,3)-x(:,3),'b',t,-sig3(:,3),'r')
xlabel('t'), ylabel('Angular Position Errors')
title('Angular Position of Satellite In Orbit')

figure (3)
plot(t,xe(:,1))
xlabel('t'), ylabel('Estimate of Position')
title('Position of Earth Orbiting Satellite')

figure(4)
plot(t,xe(:,3))
xlabel('t'), ylabel('Estimate of Angular Position' )
title('Angular Position of Satellite In Orbit')

%% Function 

function f=rungakutta(x)

% Function Routine for Satellite Orbiting Equation
f = [x(2);x(1)*((x(4))^2)-((6.6742*10^-11)*(5.98*10^24))/(x(1))^2 ;x(4);(-2*x(4)*x(2))/x(1)];

end
