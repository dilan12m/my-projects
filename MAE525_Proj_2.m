clear;
close all;
clc;

%% Problem 1

t = (0:1:7200)';
m = length(t);
x = zeros(m,1); % pre allocating memory
y = zeros(m,1); % pre allocating memory
v = 0.3927; % desired theta in unit of rads

for i = 1:m
    
tt  = t(i);
phi = 0.001745*tt; % simple integral so dod not bother doing it in matlab, for IC φ (0) = 0
psi = 0.04859*tt;  % simple integral so dod not bother doing it in matlab, for IC  ψ (0) = 0

e1 = cos(v)*cos(psi)*cos(phi)-((cos(v))^2)*sin(psi)*sin(phi)+((sin(v))^2)*sin(phi); % theta = v , s= psi, p = phi
e2 = cos(v)*cos(psi)*sin(phi)+((cos(v))^2)*sin(psi)*cos(phi)-((sin(v))^2)*cos(phi);
e3 = cos(v)*sin(v)*(sin(psi)+1);

x(i,:) = e1/(1+e3);
y(i,:) = e2/(1+e3);

end

figure(1);
plot(y,x)
xlabel('x'), ylabel('y')

%% Problem 2

t2 = (0:.1:600)';
m2 = length(t2);
x2 = zeros(m2,1); % pre allocating memory
y2 = zeros(m2,1); % pre allocating memory
w  = zeros(m2,3); % pre allocating memory
v  = 0.3927; % desired theta in unit of rads

w2  = zeros(m2,3);
phi1=zeros(m2,1);
psi1=zeros(m2,1);
v1  = ones(m2,1)*0.3927;
for j = 1:m2
    
    tt2  = t2(j);
    
    phi1(j,:) = 0.001745*tt2; % simple integral so dod not bother doing it in matlab, for IC  φ (0) = 0
  
    psi1(j,:) = 0.04859*tt2;
    
    phi = 0.001745*tt2; % simple integral so dod not bother doing it in matlab, for IC  φ (0) = 0
    psi = 0.04859*tt2;  % simple integral so dod not bother doing it in matlab, for IC  ψ (0) = 0

    e12 = cos(v)*cos(psi)*cos(phi)-((cos(v))^2)*sin(psi)*sin(phi)+((sin(v))^2)*sin(phi); % theta = v , s= psi, p = phi
    e22 = cos(v)*cos(psi)*sin(phi)+((cos(v))^2)*sin(psi)*cos(phi)-((sin(v))^2)*cos(phi);
    e32 = cos(v)*sin(v)*(sin(psi)+1);

    x2(j,:) = e12/(1+e32);
    y2(j,:) = e22/(1+e32);

    %Angular Velocity Trajectories
     w1 = inv((1/sin(v)*[cos(v)*sin(phi) -cos(v)*cos(phi)  -sin(v);...
                        -sin(v)*cos(phi) -sin(v)*sin(phi)   0     ;...
                        -sin(phi)         cos(phi)          0     ])) * [0.001745 0 0.04859 ]';
 
                    
    % Checking if Euler Angles are correct
    EulerRates = (1/sin(v)*[cos(v)*sin(phi) -cos(v)*cos(phi)  -sin(v);...
                           -sin(v)*cos(phi) -sin(v)*sin(phi)   0     ;...
                           -sin(phi)         cos(phi)          0     ])*w1;
    

    % desired angular velocity trajectories 
    w(j,:) = w1; 
    w2(j,:) =  EulerRates; 


end

figure(2);
plot(y2,x2)
xlabel('x'), ylabel('y')

figure(3);
subplot(3,1,1)
plot (t2, w(:,1))
ylabel('w1')
xlabel('t')
title ('Desired Angular Velocity')
subplot(3,1,2)
plot (t2, w(:,2))
ylabel('w2')
xlabel('t')
title ('Desired Angular Velocity')
subplot(3,1,3)
plot (t2, w(:,3))
ylabel('w3')
xlabel('t')
title ('Desired Angular Velocity')

% plots of euler rates 
figure (4)
subplot(3,1,1)
plot (t2, phi1(:,1))
title ('Euler Rates')
subplot(3,1,2)
plot (t2, v1(:,1))
title ('Euler Rates')
subplot(3,1,3)
plot (t2, psi1(:,1))
title ('Euler Rates')


%% Problem 3

% My Euler angle attitude matrix using 3-1-3 attitude sequence at the
% initial time 
deg = 22.5;% in radians here
A = [cos(0)*cos(0) - sin(0)*cos(deg)*sin(0), sin(0)*cos(0)  + cos(0)*cos(deg)*sin(0), sin(deg)*sin(0);...
    -cos(0)*sin(0) - sin(0)*cos(deg)*cos(0), -sin(0)*sin(0) + cos(0)*cos(deg)*cos(0), sin(deg)*cos(0);...
     sin(0)*sin(deg),                        -cos(0)*sin(deg),                        cos(deg)];

% initial time desired quaternion 
q1 =sqrt(.25*(1+A(1,1)-A(2,2)-A(3,3))); 
q2 =sqrt(.25*(1+A(2,2)-A(1,1)-A(3,3)));
q3 =sqrt(.25*(1+A(3,3)-A(2,2)-A(1,1)));
q4 =sqrt(.25*(1+A(1,1)+A(2,2)+A(3,3)));

q = [q1 q2 q3 q4]'


%% Problem 3 continued

Q  = zeros(m2,4);
Q(1,:)=[q];

for k = 1:m2-1

    c12    = (sin(.5*norm(w(k,:))*.1)*w(k,:)) / norm(w(k,:));
    c21    = -(c12);
    c22    = (cos(.5*sqrt(w(k,1)^2 + w(k,2)^2 + w(k,3)^2)*.1));
    c11    = cos(.5*norm(w(k,:))*.1)*eye(3) -[0 -c12(3) c12(2); c12(3) 0 -c12(1); -c12(2) c12(1) 0] ;
    
    om     = [c11 c12'; c21 c22 ]*Q(k,:)';
     
    % desired quaternion 
    Q(k+1,:) = om'; 
   
end

figure(5);
subplot(2,2,1)
plot(t2,Q(:,1));
ylabel('Desired q1')
xlabel('t')
title('Desired Quaternion')

subplot(2,2,2)
plot(t2,Q(:,2));
ylabel('Desired q2')
xlabel('t')
title('Desired Quaternion')

subplot(2,2,3)
plot(t2,Q(:,3));
ylabel('Desired q3')
xlabel('t')
title('Desired Quaternion')

subplot(2,2,4)
plot(t2,Q(:,4));
ylabel('Desired q4')
xlabel('t')
title('Desired Quaternion')

%% Problem 4

L  = zeros(m2,3);
q0 = [0 0 0 1]';
w0 = [0 0 0]';

dt = .1;
xr = zeros(m2,7);
xr(1,:)=[q0;w0]'; % not desired q and w 

for kk = 1:m2-1
    
    
    kq = 15;
    kw = 15;

    EEd = [ Q(kk,4) -Q(kk,3)  Q(kk,2);...
            Q(kk,3)  Q(kk,4) -Q(kk,1);...
           -Q(kk,2)  Q(kk,1)  Q(kk,4);...
           -Q(kk,1) -Q(kk,2) -Q(kk,3)];
    
    L(kk,:) = -kq * EEd' * xr(kk,1:4)' - kw*( xr(kk,5:7)' - w(kk,:)');

        f1=dt*rungakutta(xr(kk,:),L(kk,:));
        f2=dt*rungakutta(xr(kk,:)+0.5*f1',L(kk,:));
        f3=dt*rungakutta(xr(kk,:)+0.5*f2',L(kk,:));
        f4=dt*rungakutta(xr(kk,:)+f3',L(kk,:));
        xr(kk+1,:)=xr(kk,:)+(1/6*(f1'+2*f2'+2*f3'+f4'));
        
        

end
 
del_q = zeros(m2,3);

for k = 1:m2-1
    
     E = [ Q(k,4) -Q(k,3)  Q(k,2);...
           Q(k,3)  Q(k,4) -Q(k,1);...
          -Q(k,2)  Q(k,1)  Q(k,4);...
          -Q(k,1) -Q(k,2) -Q(k,3)];
                
     U = ((quatnorm(Q(k,:)))^-2) * E' * xr(k,1:4)';
     del_q(k+1,:) = U; 
end

figure(6);
subplot(2,2,1)
plot(t2,xr(:,1));
ylabel('q1')
xlabel('t')
title ('Quaternion')

subplot(2,2,2)
plot(t2,xr(:,2));
ylabel('q2')
xlabel('t')
title ('Quaternion')

subplot(2,2,3)
plot(t2,xr(:,3));
ylabel('q3')
xlabel('t')
title ('Quaternion')

subplot(2,2,4)
plot(t2,xr(:,4));
ylabel('q4')
xlabel('t')
title ('Quaternion')



figure(7);
subplot(2,2,1)
plot(t2,xr(:,1),'r',t2,Q(:,1));
ylabel('q1')
xlabel('t')
title ('Quaternion vs Desired Quaternion')

subplot(2,2,2)
plot(t2,xr(:,2),'r',t2,Q(:,2));
ylabel('q2')
xlabel('t')
title ('Quaternion vs Desired Quaternion')

subplot(2,2,3)
plot(t2,xr(:,3),'r',t2,Q(:,3));
ylabel('q3')
xlabel('t')
title ('Quaternion vs Desired Quaternion')

subplot(2,2,4)
plot(t2,xr(:,4),'r',t2,Q(:,4));
ylabel('q4')
xlabel('t')
title ('Quaternion vs Desired Quaternion')



figure(8);
subplot(3,1,1)
plot (t2, xr(:,5))
ylabel('w1')
xlabel('t')
title ('Angular Velocity')

subplot(3,1,2)
plot (t2, xr(:,6))
ylabel('w2')
xlabel('t')
title ('Angular Velocity')

subplot(3,1,3)
plot (t2, xr(:,7))
ylabel('w3')
xlabel('t')
title ('Angular Velocity')


figure(9);
subplot(3,1,1)
plot (t2, xr(:,5),'r', t2, w(:,1))
ylabel('w1')
xlabel('t')
title ('Angular Velocity vs Desired Angular Velocity')

subplot(3,1,2)
plot (t2, xr(:,6),'r', t2, w(:,2))
ylabel('w2')
xlabel('t')
title ('Angular Velocity vs Desired Angular Velocity')

subplot(3,1,3)
plot (t2, xr(:,7),'r',t2, w(:,3))
ylabel('w3')
xlabel('t')
title ('Angular Velocity vs Desired Angular Velocity')



figure(10);
subplot(3,1,1)
plot (t2, L(:,1))
ylabel('L1')
xlabel('t')
title ('Control Law')

subplot(3,1,2)
plot (t2, L(:,2))
ylabel('L2')
xlabel('t')
title ('Control Law')

subplot(3,1,3)
plot (t2, L(:,3))
ylabel('L3')
xlabel('t')
title ('Control Law')



figure(11);
subplot(3,1,1)
plot (t2, xr(:,5)-w(:,1))
ylabel('w1')
xlabel('t')
title ('Difference in angular velocity Desired and Not Desired ')

subplot(3,1,2)
plot (t2, xr(:,6)-w(:,2))
ylabel('w2')
xlabel('t')
title ('Difference in angular velocity Desired and Not Desired ')

subplot(3,1,3)
plot (t2,xr(:,7)- w(:,3))
ylabel('w2')
xlabel('t')
title ('Difference in angular velocity Desired and Not Desired ')

figure(12);
plot(t2,del_q(:,1),'r',t2,del_q(:,2),'g',t2,del_q(:,3),'b');
ylabel('Del q')
xlabel('t')
title ('Quaternion at each time step ')

mag = zeros(m2,1);

for i = 1:m2-1
    
    mag(i,:) = norm(L(i,:));
    
    
end

figure(13);
plot(t2,mag);
ylabel('MAG L')
xlabel('t')
title ('MAGINTUDE OF CONTROLLER L ')

function f=rungakutta(x0,L1) 
EE = [x0(4) -x0(3)  x0(2);...
      x0(3)  x0(4) -x0(1);...
     -x0(2)  x0(1)  x0(4);...
     -x0(1) -x0(2) -x0(3)];
 
WW = [x0(5) x0(6) x0(7)]';

J  = [ 399   -2.71  -1.21;...
      -2.71   377    2.14;...
      -1.21   2.14   377];  
  
WC = [0     -x0(7)  x0(6);...
      x0(7)   0    -x0(5);...
     -x0(6)  x0(5)   0];

dq = .5 * EE * WW;
dw = -inv(J)* WC * (J * WW) + inv(J)  * L1';

f=[dq ; dw];
end


