clear;
close all;
clc;

%% Problem 1 Everything left in Km units not m

% Initial Conditions 
mu = 398600.64;
R  = [5371.844151186472 -4141.062031065303 460.1400917227622]';
R0 = [3.230645742388105 3.522344029484922 -5.981911152962826]';


r = norm(R);
v = norm(R0);

a  = inv( (2/r) - ((v^2)/mu) ); % semimajor axis 
hh = cross (R,R0);
h  = norm(hh);
ee = (cross(R0,hh)/mu) - R/r;
e  = norm (ee);% eccentricity

ie = (ee/e)';
ih = (hh/h)';
ip = (cross(ih,ie));

A = [ie;ip;ih];

inclination  = acos(A(3,3)); %  inclination  

ascending_node = atan2(A(3,1),-A(3,2));  % ascending node 
periapsis  = atan2(A(1,3), A(2,3));   % periapsis
sig = (R'*R0) / sqrt(mu);

E0 = atan2( ((sig)/(sqrt(a))) , 1 - (r/a));

M0 = E0 - e*sin(E0); %% Mean anomaly 


%% Problem 2 

n = sqrt(mu/ (a^3));
M = M0 + n*(1);
E = vpa(KeplersEquation(e,M)); % solving keplers eq. for E 

rr = a*(1 - e*cos( E ));
x  = a*(cos( E ) - e);
y  = a * sqrt(1 - e^2)* sin(E );

x_dot = - ((sqrt(mu*a)) / rr ) * sin (E);
y_dot = (sqrt(mu*a*(1 - e^2)) / rr) * cos(E);

xx = [x_dot y_dot 0]'; %
yy = [x y 0]';
% Gives velocity vector 
rrr = vpa(A'*xx); % using vpa for better presion becaue ode45 was not working, but did not help
% Gives position vector  
roo = vpa(A'*yy);

%% Problem 3

t = [1:10:86400];
m = length(t);
dt =10;

xr = zeros(m,6);
xr(1,:)=[R;R0]';

for i=1:m-1
f1=dt*rungakutta(xr(i,:));
f2=dt*rungakutta(xr(i,:)+0.5*f1');
f3=dt*rungakutta(xr(i,:)+0.5*f2');
f4=dt*rungakutta(xr(i,:)+f3');
xr(i+1,:)=xr(i,:)+1/6*(f1'+2*f2'+2*f3'+f4');
end

% Matlab Solution Using ODE45 and ODE23tb
% Used options here to modoify the ode45 structure in order to fix the
% problem that was happening and it worked like a charm

options = odeset('RelTol',1e-13,'AbsTol',1e-14);
[t,xr1] = ode45(@orbitequation,t,[R;R0]', options);
[t,xr2] = ode23tb(@orbitequation,t,[R;R0]');

figure (1);
subplot(311)
plot3(xr(:,1),xr(:,2),xr(:,3))
title('Runga Kutta Method')
xlabel('x(t)')
ylabel('y(t)')
zlabel('z(t)')

subplot(312)
plot3(xr1(:,1),xr1(:,2),xr1(:,3))
title('ODE45')
xlabel('x(t)')
ylabel('y(t)')
zlabel('z(t)')

subplot(313)
plot3(xr2(:,1),xr2(:,2),xr2(:,3))
title('ODE23tb')
xlabel('x(t)')
ylabel('y(t)')
zlabel('z(t)')


%% Problem 4
rw = zeros(m,3);
rww = zeros(m,3);
ti = linspace(10,86400,8640);

for i = 1:m
M1 = M0 + n * ti(i);
E1 = (KeplersEquation(e,M1));
r1 = a * (1 - e*cos(E1));
x1 = a*(cos(E1)- e);
y1 = a*sqrt(1-(e^2))*sin(E1);
x1_dot = -((sqrt(mu*a))/(r1))*sin(E1);
y1_dot = ((sqrt(mu*a*(1-(e^2)))) / r1) * cos(E1);
rw(i,:) = (A'*[x1 y1 0]')';
rww(i,:) = (A'*[x1_dot y1_dot 0]')';
end


figure(2);
plot3(rw(:,1),rw(:,2),rw(:,3))
title('Orbit Elements Apporach')
axis equal
xlabel('x(t)')
ylabel('y(t)')
zlabel('z(t)')

%% Problem 5
rw2 = zeros(m,3);
rww2 = zeros(m,3);
sigma = (R'*R0)/(sqrt(mu));

for i = 1:m
    
j = ti(i);
ET = Etilda(r,a,sigma,mu,j);
r_t = a + (r - a)* cos(ET) + sigma*sqrt(a)*sin(ET);
F = 1 - (a/r)*(1 - cos(ET));
F_dot = - ((sqrt(mu*a))/(r_t*r))*sin(ET);
G = ti(i) + sqrt((a^3)/(mu))*(sin(ET)-ET);
G_dot = 1 - (a/r_t)*(1-cos(ET));
rw2(i,:) = F*R + G*R0;
rww2(i,:) = F_dot*R + G_dot*R0;

end


figure(3);
plot3(rw2(:,1),rw2(:,2),rw2(:,3))
title('F and G Apporach')
axis equal
xlabel('x(t)')
ylabel('y(t)')
zlabel('z(t)')


%% Function Used

function E = KeplersEquation(e,M)
f = @(E) E - e * sin(E) - M;
E = fzero(f,0);
return 
end

function E2 = Etilda(r,a,sigma,mu,j)
f = @(E2) E2 - (1 - (r/a))*sin(E2) - (sigma/sqrt(a))*(cos(E2)-1) - (sqrt((mu)/(a^3)))*j;
E2 = fzero(f,0);
return 
end

function dydt=orbitequation(t,x)
%     r=sqrt((x(1)^2)+(x(2)^2)+(x(3)^2));
      r= norm([x(1),x(2),x(3)]);
      dydt=[x(4); x(5); x(6); -(398600.64*x(1))/r^3; -(398600.64*x(2))/r^3; -(398600.64*x(3))/r^3];  
    
end

function f=rungakutta(x)

% Function Routine for
r= norm([x(1),x(2),x(3)]);
f=[x(4); x(5); x(6); -(398600.64*x(1))/r^3; -(398600.64*x(2))/r^3; -(398600.64*x(3))/r^3];

end