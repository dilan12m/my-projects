clear all; close all; clc;

% System parameter
ell = .5; % m
g = 9.81;
R = 0.01; %radius of pully at moter

% continous time state space
Ac = [0 1; 3*g/(2*ell) 0];
Bc = [3*R/(2*ell); 0];
C = [1 0];
D = 0;

%simulating continous time initial condition response
% x0 = [ell *pi/180; 0];
x0 = [1*pi/180; 0];
myConSys = ss(Ac,Bc,C,D);
tFinal = .8; % final time simulation in seconds
[yCon, t,~] = initial (myConSys, x0, tFinal);
figure
% plot(t, yCon*180/pi, 'k-' , tDis, yDis*180/pi, 'ro' ,'linewi' , 2), grid on
% xlabel('t[s]'), ylabel('\theta [deg]')

% converting our state space model to discrete time 
fs = 20;
T = 1/fs;
A = expm(Ac*T);
B = (A-eye(2))/(Ac)*Bc;
myDisSys = c2d(myConSys,T);
% A = myDisSys.A; % pullinf stuff from c2d
% B = myDisSys.B;


myDisSys2 = ss(A,B,C,D,T);
tFinal = .8; % final time simulation in seconds
[yDis, tDis,~] = initial (myDisSys, x0, tFinal);
% figure
% plot(tDis, yDis*180/pi, 'k-' , 'linewi' , 2), grid on
% xlabel('t[s]'), ylabel('\theta [deg]')


plot(t, yCon*180/pi, 'k-' , tDis, yDis*180/pi, 'ro' ,'linewi' , 2), grid on
xlabel('t[s]'), ylabel('\theta [deg]')
legend( 'continous','discrete')


%% apply an zero-order-hold to the system

x0 = [0;0];
N = 20;
t = (0:T:(N-1)*T);
u = sin(1.6*2*pi*t+ pi/4);
yCon2 = lsim(myConSys, u, t, x0,'zoh');
yDis2 = lsim(myDisSys, u, t, x0);

figure
subplot(211)
stairs(t,u,'linewi' , 2),grid on
subplot(212)
plot(t,yCon2*180/pi,'k-',...
     t,yDis2*180/pi,'ro', 'linewi' , 2),grid on


%% state feedback controller

zstar = [0.76+0.31i, 0.76-0.31i];
K = acker(A,B,zstar);
A_CL = A - B*K;
eig(A_CL);

%% simulating the closed loop system

x0 = [15*pi/180; 0];
y0 =15*pi/180 ;%initial output(initial angle of pendulum)


N =20;
t = (0:T:(N-1)*T)';

y = zeros(1,N);
x = zeros(2,N);
u = zeros(1,N);

y(1) = y0;
x(:,1) = x0;
u(1) = -K*x(:,1);

for k = 2:N
    x(:,k) = A*x(:,k-1) +B*u(k-1);
    u(k) = -K*x(:,k);
    y(k) = C*x(:,k) +D*u(k);
  
end

figure
subplot(211), stairs (t,u,'k-', 'linewi',2), grid on
title('control input motar angular velocity [degs] ');
subplot(212), plot (t,y,'k.','markersize',12), grid on
title('output pendulum angle [degs] ');


%% estimator

betastar= [0.21 + 0.51i, 0.21 - 0.51i];
L = acker(A',C',betastar)'; %acker siso and could use place works for multi input multi output
% (A-L*C) checking eig -> desired roots


%% simulating the closed loop system (w/estimator)

% x0 = [15*pi/180; 0];
% y0 =  15*pi/180 ;%initial output(initial angle of pendulum)

y0 = 0;
x0 = [0;0];

N =200; % # of time step
t = (0:T:(N-1)*T)';

y = zeros(1,N);% actual output( sensor measurement angle)
x = zeros(2,N);% true states (we do not know this)
u = zeros(1,N);

yhat = zeros(1,N);%  estimited output (angle of the pendulum)
xhat = zeros(2,N);% estimited states


y(1) = y0; % initial sensor measurement 
yhat(1) = y0;
x(:,1) = x0;
xhat(:,1) = x0;
u(1) = -K*xhat(:,1);
% cat = 5*square(1*pi*t, 0.2) - 5*square(1*pi*t-pi,0.2);
% % stairs(t,cat)
cat = zeros(1,N);
for k = 2:N

    if k == 3
         cat(k) = 500*rand(1);
        x(:,k) = A*x(:,k-1) + B*(u(k-1)+cat(k));% true system dynamics
    elseif k == 100
        cat(k) = 500*rand(1);
        x(:,k) = A*x(:,k-1) + B*(u(k-1)+cat(k));% true system dynamics
     elseif k == 150
        cat(k) = 500*rand(1);
        x(:,k) = A*x(:,k-1) + B*(u(k-1)+cat(k));% true system dynamics
    else
        x(:,k) = A*x(:,k-1) + B*u(k-1);% true system dynamic
    end
%   x(:,k) = A*x(:,k-1) + B*u(k-1);% true system dynamics
    xhat(:,k) = A*xhat(:,k-1) + B*u(k-1) - L*(yhat(k-1) - y(k-1));
    u(k) = -K*xhat(:,k);
    yhat(k) = C*xhat(:,k) + D*u(k);
    y(k) = C*x(:,k) + D*u(k) + 0.01*randn(1);
end

figure
subplot(211), stairs (t,u*180/pi,'k-', 'linewi',2), grid on
hold on
stairs(t,cat*180/pi,'r-','linewi',2)
title('control input motar angular velocity [degs] ');
legend('control input', 'cat' )
subplot(212), plot (t,y*180/pi,'k.',t,yhat*180/pi,'g.','markersize',12), grid on
title('output pendulum angle [degs] ');
legend('y measured angle', 'yhat estimited measurement')

%% first calculate position of cart at each time step

x0 = zeros(1,N);
xTop= zeros(1,N);
yTop= zeros(1,N);
xTop(1) = x0(1) - ell*sin(x(1,1));
yTop(1) = ell*cos(x(1,1));
x0(1) = 0;

for k=2:N
    x0(k) = x0(k-1) + R*T*u(k-1);
    xTop(k) = x0(k) - ell*sin(x(1,k));
    yTop(k) =ell*cos(x(1,k));
    
end

figure 
plot (t, x0, 'linewi', 2), grid minor 
xlabel('t [s]')
legend('position of the cart [m]')

%% actual animation

for k=1:N
    figure(100)
    plot(x0(k),0,'ks', 'markersize', 25, 'linewi', 4);
    hold on;
    plot([x0(k), xTop(k)],[0,yTop(k)],'k-', 'linewi', 4);
    hold off;
    axis([-1 1 -1 1])
    axis equal
    drawnow 
    pause(T);
end



























