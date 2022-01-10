clear;
close all;
clc;
%spring mass position control using state feedback
T = .0126; %sampling period
%SS model in modal coordinites
A = [.9798 .1232; -.1232 .9798];
B = [.0127; .0375 ];
C = [.0306 , -.0083];
D = [0];

N = [C 0; A-eye(2) B]\[1; zeros(2,1)];
Nx = N(1:2);
Nu = N(3);

% r = 0.5;
% xr = Nx *r
% uss = Nu*r

%% lets run a controller

% settling time of 1s and overshoot of 10%
zstar = [0.9487 + 0.0652i, 0.9487 - 0.0652i];
K = acker(A, B, zstar);

% to check if get back our zstar
eig(A - B*K);

t = (0:T:1.5)';
N = length(t);
x = zeros(2,N);
y = zeros(1,N);
u = zeros(1,N);
r = .5*ones(1,N);
xr = Nx*r; % give all referance state
uss = Nu*r; % what we want to get to

u(1) = -K*(x(:,1) - xr(:,1));
 
for k = 2:N
    x(:,k) = A*x(:,k-1)+ B*u(k-1);
    y(k) = C*x(:,k) + D*u(k);
    u(k) = -K*(x(:,k) - xr(:,k));
end
    
figure
plot(t,r,'r-', t,y,'k.', 'markersize', 14)
xlabel('t[s]'),legend('referance position','output position')
 

%% Augumented model(adding a integral state)
    
A_aug = [A, zeros(2,1);...
        -C, 1];
    
B_aug = [B;-D];


% getting state feed matrix
zstar_aug = [zstar .8*abs(zstar(1))]';
K_aug = acker(A_aug,B_aug,zstar_aug);
eig(A_aug - B_aug*K_aug);

% define estimator gain matrix
betastar = [0.8225+ 0.2501i; 0.8225-0.2501i];
L = acker(A',C',betastar)';
eig(A - L*C);

% redefine (reset) Variable
x  = zeros(2,N);
xhat  = zeros(2,N); % estimated states
xI = zeros(1,N); % integral state
y  = zeros(1,N);
yhat  = zeros(1,N);% estimated outputs
u  = zeros(1,N);
e  = zeros(1,N);

r = .5*ones(1,N);
xr = Nx*r; % give all referance state
uss = Nu*r; % what we want to get to

u(1) = -K_aug*([xhat(:,1); xI(1)]- [xr(:,1);0]);
e(1) = r(1) - y(1);
for k = 2:N
    x(:,k) = A*x(:,k-1)+ B*u(k-1);
    xhat(:,k) = A*xhat(:,k-1)+ B*u(k-1) - L*(yhat(k-1)-y(k-1));
    xI(k)  = xI(k-1)  + e(k-1);
    u(k)   = -K_aug*([xhat(:,k);xI(k)] - [xr(:,k);0]);
    y(k)   = C*x(:,k) + D*u(k);%+0.01*randn(1);
    yhat(k)= C*xhat(:,k) + D*u(k);
    e(k)   = r(k)-y(k);
    
end
    
figure(2)
plot(t,r,'r-', t,y,'k.', 'markersize', 14)
xlabel('t[s]'),legend('referance position','output position')




    
% 
% K = acker(A , B, [0.9878 + 0.0059i , 0.9878 - 0.0059i])
% L = acker(A' , C', [ 0.9440 + 0.0875i ,  0.9440 - 0.0875i])'

    
    