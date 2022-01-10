clear;
close all;
clc;
%% USER DEFINED AREA

initial_axes = 131; 
phi = 5*pi/6; theta = pi/6; psi = pi;
qq = quaternion([rad2deg(psi),rad2deg(theta),rad2deg(phi)],'eulerd','XZX','point');


%%  Euler angle rotation

% Define initial frame
x_0 = [1;0;0]; 
y_0 = [0;1;0];
z_0 = [0;0;1];

o_v= [x_0'; y_0'; z_0']; 

% plotting axes
figure (1)
xlabel('x')
ylabel('y')
zlabel('z')
plot3(x_0(1),y_0(2),z_0(3),'b');
hold on
plot3([0;o_v(1,1)],[0;o_v(1,2)],[0;o_v(1,3)],'b')
plot3([0;o_v(2,1)],[0;o_v(2,2)],[0;o_v(2,3)],'b')
plot3([0;o_v(3,1)],[0;o_v(3,2)],[0;o_v(3,3)],'b') 

hold on

x_1 = DCM(initial_axes,phi,0,0)*x_0; 
y_1 = DCM(initial_axes,phi,0,0)*y_0;
z_1 = DCM(initial_axes,phi,0,0)*z_0;
o_1 = [x_1'; y_1'; z_1'];
plot3([0;o_1(1,1)],[0;o_1(1,2)],[0;o_1(1,3)],'r')
plot3([0;o_1(2,1)],[0;o_1(2,2)],[0;o_1(2,3)],'r')
plot3([0;o_1(3,1)],[0;o_1(3,2)],[0;o_1(3,3)],'r')
hold on 

x_2 = DCM(initial_axes,0,theta,0)*x_1;
y_2 = DCM(initial_axes,0,theta,0)*y_1;
z_2 = DCM(initial_axes,0,theta,0)*z_1;
o_2 = [x_2'; y_2'; z_2'];
plot3([0;o_2(1,1)],[0;o_2(1,2)],[0;o_2(1,3)],'g')
plot3([0;o_2(2,1)],[0;o_2(2,2)],[0;o_2(2,3)],'g')
plot3([0;o_2(3,1)],[0;o_2(3,2)],[0;o_2(3,3)],'g')

x_3 = DCM(initial_axes,0,0,psi)*x_2;
y_3 = DCM(initial_axes,0,0,psi)*y_2;
z_3 = DCM(initial_axes,0,0,psi)*z_2;
o_3 = [x_3'; y_3'; z_3']
plot3([0;o_3(1,1)],[0;o_3(1,2)],[0;o_3(1,3)],'c','linewi',2)
plot3([0;o_3(2,1)],[0;o_3(2,2)],[0;o_3(2,3)],'c','linewi',2)
plot3([0;o_3(3,1)],[0;o_3(3,2)],[0;o_3(3,3)],'c','linewi',2)
plot3(o_3(1,1),o_3(1,2),o_3(1,3),'bd')
plot3(o_3(2,1),o_3(2,2),o_3(2,3),'rd')
plot3(o_3(3,1),o_3(3,2),o_3(3,3),'gd')
axis([-1 1 -1 1 -1 1])
xlabel('x')
ylabel('y')
zlabel('z')
grid on
title('Eular Angles')


%% Quaternions
a = [1,0,0];
b = [0,1,0];
c = [0,0,1];
% qq = quaternion([rad2deg(psi),rad2deg(theta),rad2deg(phi)],'eulerd','XZX','point');
rP = rotmat(qq,'point')


figure (2);
plot3(a(1),a(2),a(3),'b');

hold on
grid on
axis([-1 1 -1 1 -1 1])
xlabel('x')
ylabel('y')
zlabel('z')

plot3(b(1),b(2),b(3),'r');
plot3(c(1),c(2),c(3),'c');
plot3(rP(1,1),rP(1,2),rP(1,3),'bd')
plot3(rP(2,1),rP(2,2),rP(2,3),'rd')
plot3(rP(3,1),rP(3,2),rP(3,3),'gd')

plot3([0;rP(1,1)],[0;rP(1,2)],[0;rP(1,3)],'c','linewi',2)
plot3([0;rP(2,1)],[0;rP(2,2)],[0;rP(2,3)],'c','linewi',2)
plot3([0;rP(3,1)],[0;rP(3,2)],[0;rP(3,3)],'c','linewi',2)

plot3([0;a(1)],[0;a(2)],[0;a(3)],'b')
plot3([0;b(1)],[0;b(2)],[0;b(3)],'b')
plot3([0;c(1)],[0;c(2)],[0;c(3)],'b')
title('Quaternion')




%% Exponential Coordinates for Rotation 

veda = ( acos((trace(o_3)-1)/2));
omega = (1/(2*sin(veda))) * [o_3(3,2)-o_3(2,3); o_3(1,3)- o_3(3,1);o_3(2,1)- o_3(1,2)];

X = [  0         -omega(3)   omega(2) ;
       omega(3)   0         -omega(1) ;
      -omega(2)   omega(1)   0       ];

omega2 = omega*omega' -  ((norm (omega))^2)*eye(3);

if norm (omega) == 1
    expo1 = eye(3) + X * sin(veda) + omega2 * (1-cos(veda))
    
else 
    expo1 = eye(3) + (X/norm (omega)) * sin(norm(omega)*veda) + (omega2/(norm(omega))^2) * (1-cos(norm(omega)*veda))
    
end

figure (3);
plot3(a(1),a(2),a(3),'b');

hold on
grid on
axis([-1 1 -1 1 -1 1])
xlabel('x')
ylabel('y')
zlabel('z')

plot3(b(1),b(2),b(3),'r');
plot3(c(1),c(2),c(3),'c');
plot3(expo1(1,1),expo1(1,2),expo1(1,3),'bd')
plot3(expo1(2,1),expo1(2,2),expo1(2,3),'rd')
plot3(expo1(3,1),expo1(3,2),expo1(3,3),'gd')

plot3([0;expo1(1,1)],[0;expo1(1,2)],[0;expo1(1,3)],'c','linewi',2)
plot3([0;expo1(2,1)],[0;expo1(2,2)],[0;expo1(2,3)],'c','linewi',2)
plot3([0;expo1(3,1)],[0;expo1(3,2)],[0;expo1(3,3)],'c','linewi',2)

plot3([0;a(1)],[0;a(2)],[0;a(3)],'b')
plot3([0;b(1)],[0;b(2)],[0;b(3)],'b')
plot3([0;c(1)],[0;c(2)],[0;c(3)],'b')
title('Exponential Coordinates for Rotation')


%%





%Attitude matrix
function A = DCM(euler_rotation, psi, theta, phi)

    if euler_rotation==121
        A = [cos(theta)          sin(phi)*sin(theta)                             -cos(phi)*sin(theta);
             sin(theta)*sin(psi) cos(phi)*cos(psi)-sin(phi)*cos(theta)*sin(psi)  sin(phi)*cos(psi)+cos(phi)*cos(theta)*sin(psi);
             sin(theta)*cos(psi) -cos(phi)*sin(psi)-sin(phi)*cos(theta)*cos(psi) -sin(phi)*sin(psi)+cos(phi)*cos(theta)*cos(psi)];
    
     
    elseif euler_rotation==232
        A = [-sin(phi)*sin(psi)+cos(phi)*cos(theta)*cos(psi)  sin(theta)*cos(psi)  -cos(phi)*sin(psi)-sin(phi)*cos(theta)*cos(psi);
             -cos(phi)*sin(theta)                             cos(theta)            sin(phi)*sin(theta);
             sin(phi)*cos(psi)+cos(phi)*cos(theta)*sin(psi)   sin(theta)*sin(psi)   cos(phi)*cos(psi)-sin(phi)*cos(theta)*sin(psi) ];
    
    elseif euler_rotation==313
        A = [cos(phi)*cos(psi)-sin(phi)*cos(theta)*sin(psi)   sin(phi)*cos(psi)+cos(phi)*cos(theta)*sin(psi)  sin(theta)*sin(psi);
             -cos(phi)*sin(psi)-sin(phi)*cos(theta)*cos(psi)  -sin(phi)*sin(psi)+cos(phi)*cos(theta)*cos(psi) sin(theta)*cos(psi);
             sin(phi)*sin(theta)                              -cos(phi)*sin(theta)                            cos(theta)];
  
    elseif euler_rotation==131
         A = [cos(theta)          cos(phi)*sin(theta)                              sin(phi)*sin(theta);
             -sin(theta)*cos(psi) -sin(phi)*sin(psi)+cos(phi)*cos(theta)*cos(psi)  cos(phi)*sin(psi)+sin(phi)*cos(theta)*cos(psi);
             sin(theta)*sin(psi)  -sin(phi)*cos(psi)-cos(phi)*cos(theta)*sin(psi)  cos(phi)*cos(psi)-sin(phi)*cos(theta)*sin(psi)];
        
    elseif euler_rotation==212
        A = [cos(phi)*cos(psi)-sin(phi)*cos(theta)*sin(psi) sin(theta)*sin(psi)  -sin(phi)*cos(psi)-cos(phi)*cos(theta)*sin(psi);
             sin(phi)*sin(theta)                            cos(theta)           cos(phi)*sin(theta);
             cos(phi)*sin(psi)+sin(phi)*cos(theta)*cos(psi) -sin(theta)*cos(psi) -sin(phi)*sin(psi)+cos(phi)*cos(theta)*cos(psi)];

    elseif euler_rotation==323
         A = [-sin(phi)*sin(psi)+cos(phi)*cos(theta)*cos(psi) cos(phi)*sin(psi)+sin(phi)*cos(theta)*cos(psi)  -sin(theta)*cos(psi);
              -sin(phi)*cos(psi)-cos(phi)*cos(theta)*sin(psi) cos(phi)*cos(psi)-sin(phi)*cos(theta)*sin(psi)   sin(theta)*sin(psi);
              cos(phi)*sin(theta)                             sin(phi)*sin(theta)                              cos(theta)];
     
    elseif euler_rotation==123
        A = [cos(theta)*cos(psi)  cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi) sin(phi)*sin(psi)-cos(phi)*sin(theta)*cos(psi);
            -cos(theta)*sin(psi)  cos(phi)*cos(psi)-sin(phi)*sin(theta)*sin(psi) sin(phi)*cos(psi)+cos(phi)*sin(theta)*cos(psi);
            sin(theta)            -sin(phi)*cos(theta)                           cos(phi)*cos(theta)];
       
    elseif euler_rotation==231
        A = [cos(phi)*cos(theta)                                sin(theta)           -sin(phi)*cos(theta) ;
            sin(phi)*sin(psi)-cos(phi)*sin(theta)*cos(psi)      cos(phi)*cos(psi)    cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi)
            sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi)      -cos(phi)*sin(psi)   cos(phi)*cos(psi)-sin(phi)*sin(theta)*sin(psi)];
  
    elseif euler_rotation==312
        A = [cos(phi)*cos(psi)-sin(phi)*sin(theta)*sin(psi)   sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi)   -cos(theta)*sin(psi);
            -sin(phi)*cos(theta)                              cos(phi)*cos(theta)                               sin(theta);
            cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi)    sin(phi)*sin(psi)-cos(phi)*sin(theta)*cos(psi)    cos(theta)*cos(psi)];
        
    elseif euler_rotation==132
        A =[ cos(theta)*cos(psi)  sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi) -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi);
            -sin(theta)           cos(phi)*cos(theta)                             sin(phi)* cos(theta);
            cos(theta)*sin(psi)   -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi) cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi)];
        
    elseif euler_rotation==213
        A = [cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi)  cos(theta)*sin(psi) -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi);
             -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi) cos(theta)*cos(psi)  sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi);
             sin(phi)*cos(theta)                             -sin(theta)          cos(phi)*cos(theta)];
       
    elseif euler_rotation==321
        A = [cos(phi)*cos(theta)                                sin(phi)*cos(theta)                                 -sin(theta) ;
             -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi)    cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi)      cos(theta)*sin(psi);
              sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi)    -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi)     cos(theta)*cos(psi)];
    
    else
        error(' Invalid input, Try again!!');    
    end
end

% 
% %Euler angle to quaternion conversions
% function q = quat(rotation_axes,psi, theta, phi)
% 
%     c1=cos(psi/2); c2=cos(theta/2); c3=cos(phi/2);
%     s1=sin(psi/2); s2=sin(theta/2); s3=sin(phi/2);
%     
%     c13 = cos((psi+phi)/2);  s13 = sin((psi+phi)./2);
%     
%     c1_3 = cos((psi-phi)/2);  s1_3 = sin((psi-phi)/2);
%     c3_1 = cos((phi-psi)/2);  s3_1 = sin((phi-psi)/2);
% 
%     if rotation_axes == 121
%         q =[c2*s13,s2*c1_3,s2*s1_3,c2*c13];
%         
%     elseif rotation_axes == 232
%         q =[s2*s1_3,c2*s13,s2*c1_3,c2*c13];
%         
%     elseif rotation_axes == 313
%         q =[s2*c1_3,s2*s1_3,c2*s13,c2*c13];
%         
%     elseif rotation_axes == 131
%         q =[c2*s13,s2*s3_1,s2*c3_1,c2*c13];
%         
%     elseif rotation_axes == 212
%         q =[s2*c3_1,c2*s13,s2*s3_1,c2*c13];
%         
%     elseif rotation_axes == 323
%         q =[s2*s3_1,s2*c3_1,c2*s13,c2*c13];
%         
%     elseif rotation_axes == 123
%         q =[s1*c2*c3+c1*s2*s3,c1*s2*c3-s1*c2*s3,c1*c2*s3+s1*s2*c3,c1*c2*c3-s1*s2*s3];
%         
%     elseif rotation_axes == 231
%         q =[c1*c2*s3+s1*s2*c3,s1*c2*c3+c1*s2*s3,c1*s2*c3-s1*c2*s3,c1*c2*c3-s1*s2*s3];
%         
%     elseif rotation_axes == 312
%         q =[c1*s2*c3-s1*c2*s3,c1*c2*s3+s1*s2*c3,s1*c2*c3+c1*s2*s3,c1*c2*c3-s1*s2*s3];
%         
%     elseif rotation_axes == 132
%         q =[s1*c2*c3-c1*s2*s3,c1*c2*s3-s1*s2*c3,c1*s2*c3+s1*c2*s3,c1*c2*c3+s1*s2*s3];
%         
%     elseif rotation_axes == 213
%         q =[c1*s2*c3+s1*c2*s3,s1*c2*c3-c1*s2*s3,c1*c2*s3-s1*s2*c3,c1*c2*c3+s1*s2*s3];
%         
%     elseif rotation_axes == 321
%         q =[c1*c2*s3-s1*s2*c3,c1*s2*c3+s1*c2*s3,s1*c2*c3-c1*s2*s3,c1*c2*c3+s1*s2*s3];
%     else
%         error(' Invalid input, Try again!!');    
%     end
% end
% 
% 
% 

% psi = pi/6;
% x_3 = DCM(initial_axes,0,0,psi)*x_2;
% y_3 = DCM(initial_axes,0,0,psi)*y_2;
% z_3 = DCM(initial_axes,0,0,psi)*z_2;
% o_3 = [x_3'; y_3'; z_3']
% % quiver3(zeros(3,1), zeros(3,1), zeros(3,1), o_2(:,1), o_2(:,2), o_2(:,3),'c')
% plot3([0;o_3(1,1)],[0;o_3(1,2)],[0;o_3(1,3)],'c')
% plot3([0;o_3(2,1)],[0;o_3(2,2)],[0;o_3(2,3)],'c')
% plot3([0;o_3(3,1)],[0;o_3(3,2)],[0;o_3(3,3)],'c')





% % q = quat(initial_axes, psi, theta, phi)
% % % q = quat(initial_axes,0,pi/6, pi/6)
% % quat_form = quaternion(q(1),q(2),q(3),q(4));
% % % quat_to_rotation_matrix = rotmat(quat_form, 'frame')
% % %rP = rotmat(quat_form, 'frame')
% % a = [1,0,0];
% % b = [0,1,0];
% % c = [0,0,1];
% % 
% % % q1 = sqrt(.25*(1+rP(1,1)-rP(2,2)-rP(3,3)));
% % % q2 = sqrt(.25*(1+rP(2,2)-rP(1,1)-rP(3,3)));
% % % q3 = sqrt(.25*(1+rP(3,3)-rP(2,2)-rP(1,1)));
% % % q4 = sqrt(.25*(1+rP(1,1)+rP(2,2)+rP(3,3)));
% % 
% % 
% % qq = quaternion([0,30,150],'eulerd','XZX','point')
% % % rP1 = rotateframe(quat_form,[a;b;c])
% % 
% % rP = rotmat(qq,'point')









