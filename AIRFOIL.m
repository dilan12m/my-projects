% 1. Generating NACA Airfoil From Equations 
% 2. Generating Reynolds Number 
% 3. Simulating Laminar Flow Around Airfoil
clear ; close all; clc; 


%% 1. Generating NACA Airfoil From Equations %%
% Referance for Mathematical Equations http://www.aerospaceweb.org/question/airfoils/q0041.shtml


% 1. Pick the value of x from 0 to max chord. 

% Chord
x = 1; 

t = [0:.01:x];
mm = length(t);

% Specifiy the 4 digit NACA foil number
NACA_Number = [2 4 15]; % enter the 4 digit NACA foil number

% Max Camber
m = NACA_Number(1)/100;

% Position Of Max Camber
p = NACA_Number(2)/10;

% Thickness Of Airfoil
tt = NACA_Number(3)/100;

% Storing Values 
y_ca    = zeros(1,mm);
x_upper = zeros(1,mm);
y_upper = zeros(1,mm);
x_lower = zeros(1,mm);
y_lower = zeros(1,mm);


if m == 0 && p==0
    % Case When The Airfoil is Symmetric
    for i = 1:mm 
      y_th = (tt/0.2)*(0.2969*sqrt(t(i)) - 0.126*t(i) - 0.3516*t(i)^2 + 0.2843*t(i)^3 - 0.1015*t(i)^4);
      x_upper(i) = t(i) ;
      y_upper(i) = y_th ;
      x_lower(i) =   x_upper(i);
      y_lower(i) = - y_upper(i);
    end
    
else
   % Case When The Airfoil is Not Symmetric
 
    for i = 1:mm
        
        % 2. Calculating the Mean Camber Line  
        if t(i)<=p
        
           y_ca(i)  = (m/p^2)*(2*p*t(i) - t(i)^2);
           theta = atan((m/p^2)*(2*p - 2*t(i)));
        else 
            
           y_ca(i) = (m/(1-p)^2)*((1-2*p)+2*p*t(i)-t(i)^2);
           theta = atan((m/(1-p)^2)*((2*p-2*t(i))));
        end
        
        % 3. Calculating the Thickness Distribution
        y_th = (tt/0.2)*(0.2969*sqrt(t(i)) - 0.126*t(i) - 0.3516*t(i)^2 + 0.2843*t(i)^3 - 0.1015*t(i)^4);
        
        % 4. Final Coordinates of upper surface

        x_upper(i) = t(i)    - y_th * sin (theta);
        y_upper(i) = y_ca(i) + y_th * cos (theta);
        x_lower(i) = t(i)    + y_th * sin (theta);
        y_lower(i) = y_ca(i) - y_th * cos (theta);

           
    end
end

figure(1);
plot(x_upper,y_upper,x_lower,y_lower,'MarkerFaceColor','c')
hold on
plot(t,y_ca,'*')  
axis equal
grid on




%% 2. Generating Reynolds Number for Airfoil %%
% Referance: https://en.wikipedia.org/wiki/Reynolds_number, http://airfoiltools.com/calculator/reynoldsnumber
% The Reynolds number is a dimensionless value that measures the ratio of inertial forces to viscous forces and 
% descibes the degree of laminar or turbulent flow. Systems that operate at the same Reynolds number will have 
% the same flow characteristics even if the fluid, speed and characteristic lengths vary.


% Velocity of the Fluid
v = 30; % in m/s

% Chord Width
l = 0.2; % in m

% Kinematic Viscosity of the Fluid ( Fluid Property, Different for diffrent Fluids) 
V = 1.4207*10^-5;

reynolds_number = (v * l) / V;






















