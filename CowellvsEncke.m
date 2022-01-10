%Determination of spacecraft’s heliocentric position and velocity 
%after time ’t’ disturbed by a planet of gravity constant ’mu3’ and 
%’orb’ orbit around the sun, by Encke’s method
% clear;
% close all;
% clc;

global orb;
global mu; 
mu=1.32712440018e11; %sun’s gravity constant

% Earth
global mu3; 
mu3=398600.4; %earth’s gravity constant
orb = [149597870;0.01667;0;0;0;-100*24*3600];%earth’s orbital elements 

% % Jupiter
% global mu3
% mu3 = 126686534;
% orb = [778547200,0.0489,deg2rad(1.305),deg2rad(274.3), deg2rad(100.4),-100*9.9259*3600]'; 

% % Venus
% global mu3
% mu3 = 324859;
% orb = [108208927.03,0.00678,deg2rad(3.3947),deg2rad(54.9), deg2rad(76.7),-100*2802*3600]'; 

% % Uranus 
% mu3 = 5793939;
% orb = [2876679082,0.04726,deg2rad(0.773),deg2rad(96.9), deg2rad(74.02),-100*17.24*3600]'; 

t=0;
i=1;
R0=1e8*[-0.27;1.475;0.001]; %spacecraft initial position
V0=[-33;-10;1]; %spacecraft’s initial velocity
R=R0;V=V0;
dt=24*3600;
tf=100*dt;

alpha = zeros(3,1);
beta  = zeros(3,1);
rbb = zeros(101,3);
vbb = zeros(101,3);


while t<=tf

    [Rb,Vb]=trajE(mu,t,R,V,t+dt);%traj of spacecraft where the sun is the center
    rbb(i,:) = Rb';
    vbb(i,:) = Vb';
    rb=norm(Rb);
    r=norm(R);
    [R3,V3]=orbit(mu,orb,t);
    ad=disturb(mu3,R,R3);
    beta=mu*dt*(1-(rb/r)^3)/rb^3+ad*dt;
    alpha=beta*dt;
    R=Rb+alpha;
    V=Vb+beta;
    Rs(:,i)=R;
    Vs(:,i)=V;
    Ts(i,1)=t;
    t=t+dt;
    i=i+1;
end

% % 
% for i = 
%     
% end



%Determination of spacecraft’s heliocentric position and velocity
%after time ’t’ disturbed by a planet of gravity constant ’mu3’ and
%’orb’ orbit around the sun, by direct numerical integration (Cowell’s method) 
options=odeset('RelTol',1e-8);
[T,X]=ode45(@perturb,[0 100*24*3600],[R0;V0],options);

 figure
% plot3(X(:,1),X(:,2),X(:,3))
% hold on
% plot3(Rs(1,:),Rs(2,:),Rs(3,:))
plot(Ts,Rs(1,:),Ts,Rs(2,:),Ts,Rs(3,:))
hold on 
plot(T,X(:,1),T,X(:,2),T,X(:,3))
% hold on
%  plot(Ts,rbb(:,1),Ts,rbb(:,2),Ts,rbb(:,3))
% plot3(X(:,1),X(:,2),X(:,3))%%%%%%%%%ignore maybe

% % % function xdot = perturb(t,x) 
% % % %(c) 2006 Ashish Tewari
% % % %     global mu;
% % %     global mu3;
% % %     global orb;
% % %     %program for the perturbed equations of motion
% % %     R=x(1:3,1); %position of s/c relative to primary
% % %     r=norm(R);
% % %     xdot(1:3,1)=x(4:6,1);
% % %     [R3,V3]=orbit(mu3,orb,t); %position and velocity of third body 
% % %     ad=disturb(mu3,R,R3); %disturbing acceleration 
% % %     xdot(4:6,1)=-mu3*R/r^3+ad;
% % % end

function xdot = perturb(t,x) 
%(c) 2006 Ashish Tewari
    global mu;
    global mu3;
    global orb;
    %program for the perturbed equations of motion
    R=x(1:3,1); %position of s/c relative to primary
    r=norm(R);
    xdot(1:3,1)=x(4:6,1);
    [R3,V3]=orbit(mu,orb,t); %position and velocity of third body 
    ad=disturb(mu3,R,R3); %disturbing acceleration 
    xdot(4:6,1)=-mu*R/r^3+ad;
end


function a=disturb(mu3,R,R3)
    %Program for calculating the disturbance acceleration ’a’ caused by a third 
    %body on a two-body orbit
    %mu3: gravitational constant of the disturbing body ’m3’
    %R: position of mass ’m2’ relative to primary mass ’m1’
    %R3: position of
    %(c) 2006 Ashish 
    r=norm(R);
    r3=norm(R3);
    R23=R-R3; 
    r23=norm(R23); 
    fx=(r23/r3)^3-1;
    a=-mu3*(R+fx*R3)/r23^3;

end


function [r,v] = orbit(mu,orb,t)
    %program for position and velocity of a body in ’orb’ elliptical orbit 
    %Elements of ’orb’: 1x1: a; 2x1:e; 3x1:i; 4x1:w; 5x1: Om; 6x1: tau 
    %(c) 2006 Ashish Tewari
    a=orb(1);
    e=orb(2);
    i=orb(3);
    w=orb(4);
    Om=orb(5);
    tau=orb(6); 
    n=sqrt(mu/a^3);
    M=-n*tau;
    E=kepler(e,M);
    r0=a*(1-e*cos(E));
    R0=a*[cos(E)-e;sqrt(1-e^2)*sin(E);0]; 
    V0=sqrt(mu*a)*[-sin(E);sqrt(1-e^2)*cos(E);0]/r0; 
    [r,v]=trajE(mu,0,R0,V0,t); %%[R,V]
%     if abs(i)>=1e-6 %%%%%%%%%%%%%%%
%         C=rotation(i,Om,w);
%     else
%         C=eye(3);
%     end
%     r=C*R;
%     v=C*V;%%%%%%%%%%%%%%

end

% 
% function E = kepler(e,M)
% f = @(E) E - e * sin(E) - M;
% E = fzero(f,0);
% return 
% end



function [R,V]=trajE(mu,t0,R0,V0,t) 
eps=1e-10;
r0=norm(R0);%
v0=norm(V0); %
alpha=dot(R0,V0);
H=cross(R0,V0);
h=norm(H);
p=h*h/mu; 
ecv=cross(V0,H)/mu-R0/r0; %
e=norm(ecv);%
ecth0=p/r0-1; 
esth0=norm(cross(ecv,R0))/r0; 
if abs(ecth0)>=eps
    th0=atan(esth0/ecth0);
if ecth0<0
    if esth0>=0
     th0=th0+pi;
   end
elseif esth0<0
    th0=th0+2*pi;
end
elseif esth0>=0 
    th0=pi/2;
else
    th0=3*pi/2;
end
ainv=-(v0*v0)/mu+2/r0;
a=1/ainv;
n=sqrt(mu/a^3); 
E0=2*atan(sqrt((1-e)/(1+e))*tan(0.5*th0)); 
tau=t0+(-E0+e*sin(E0))/n;
M=n*(t-tau);
E=kepler(e,M);
r=a*(1-e*cos(E));
f=1+a*(cos(E-E0)-1)/r0;%
g=a*alpha*(1-cos(E-E0))/mu+r0*sqrt(a/mu)*sin(E-E0); %
fd=-sqrt(mu*a)*(sin(E-E0))/(r*r0);%
gd=1+a*(cos(E-E0)-1)/r;%
R=f*R0+g*V0;%
V=fd*R0+gd*V0;%
end




function E=kepler(e,M) 
%(c) 2006 Ashish Tewari 
E=M+e*sin(M); 
fE=E-e*sin(E)-M; 
fpE=1-e*cos(E); 
dE=-fE/fpE;
E=E+dE;
eps=1e-10; %tolerance 
i=0;
 while abs(fE)>eps 
        fE=E-e*sin(E)-M;
        fpE=1-e*cos(E);
        dE=-fE/fpE;
    E=E+dE;
    i=i+1; %iteration number 
 end
end





