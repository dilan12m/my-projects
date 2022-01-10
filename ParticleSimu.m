% For the particle simulation take into account the speed of the blue particle by taking account the distance it travels and the time it takes use frames per second….. figure it out 


clear ; close all; clc; 
movielength = 30; fps = 15; % length (seconds) and fps of the output viedo
moviename = 'Ribbon';

N=2116;  % number of particles
D=.6;    % diameter of particles
DD = .25;
M=.1;    % mass of particles

Lx=20*D; % size of box
Ly=20*D;

[x y]=ndgrid(D/2:DD:Lx-D/2,D/2:DD:Ly-D/2); % place particles on grid
ii=randperm(numel(x),N);                   % N random position on grid
x=x(ii);                                   % set x-positions
y=y(ii);                                   % set x-positions

dt=.01;         % Time Step
t =[0:dt:2*pi]; % Time Vector
m = length(t);  

% Size of the Blue Particle 
p1size = 9;

%setup figure, titles, labels, ...
anim = figure;
title('Particle Simulation', 'fontsize', 12); xlabel('x'); ylabel('y');
ax = gca; ax.FontSize = 12;

% plot commands will add lines to an existing figure
ax.NextPlot = 'add';

xlim([0, Lx]);
ylim([0, Ly]);
ax.DataAspectRatio = [1 1 1];

% Initial Velocity 
vx=0; 
vy=0;

% Initial Accelertation 
ax_old=0*x;    
ay_old=0*y;

xs =zeros(m,N);   % x-position
ys =zeros(m,N);   % y-position
vxs=zeros(m,N);   % x-velocity
vys=zeros(m,N);   % y-velocity

Fx=zeros(1,N);    % Force in x
Fy=zeros(1,N);    % Force in y

for i = 1:m
% The path the blue particle is going to follow with x being t and y being p1
      p1 = 5*sin(2*t(i))+5;
      
% Boundry conditions for the Wall, for repeling the balls off the wall   
        for k =1:N
            if x(k) <= 0
              vx(k) = -vx(k);
            end
            if x(k) >= Lx
              vx(k) = -vx(k);
            end
            if y(k) >= Ly
              vy(k) = -vy(k);
            end
            if  y(k) <= 0
              vy(k) = -vy(k);
            end
        end

% Boundry conditions for the particles, for repeling the particles off of
% one another  
        for ii=1:N 
            for jj=ii+1:N    
                dy=y(ii)-y(jj); % checking the y distance 
                dx=x(ii)-x(jj); % checking the x distance 
                if sqrt(dx^2+dy^2) <= .11 %2*.01 %
                    vx(ii) = - vx(ii);
                    vy(ii) = - vy(ii);
                    vx(jj) = - vx(jj);
                    vy(jj) = - vy(jj);
                end   
            end   
        end

    for kk=1:N                     
          dyp = p1 - y(kk);
          dxp = t(i) - x(kk);
              if sqrt(dxp^2+dyp^2) <= .1
                vx(kk) = - vx(kk);
                vy(kk) = - vy(kk);   
              end
    end   
  



% Update the position from the previous information
      x=x+vx*dt+ax_old*dt^2/2;  
      y=y+vy*dt+ay_old*dt^2/2;
% For storing the updated information later used for plotting     
      xs(i,:)=x;                         
      ys(i,:)=y;

% Boundry conditions for the particle to repel away from the blue particle
% not currently working%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for kk=1:N                     
          dyp = p1 - y(kk);
          dxp = t(i) - x(kk);
%             if sqrt(dxp^2+dyp^2) <= .1
%                 vx(kk) = - vx(kk);
%                 vy(kk) = - vy(kk);
            if (-.4 <= dyp) && (dyp <=1)
                  if (-.4 <= dxp) && (dxp <=.7) 
                        dnm=sqrt(dxp.^2+dyp.^2);           
                        F=-(D/dnm);         
                        Fx(kk)=Fx(kk)+F.*dxp;                   
                        Fy(kk)=Fy(kk)+F.*dyp;  
% %                         for nn=1:N 
% %                           for mm=nn+1:N  
% %                             Fx(nn)= F.*dx;     
% %                             Fx(mm)= - F.*dx;     
% %                             Fy(nn)= F.*dy;    
% %                             Fy(mm)=- F.*dy; 
% %                           end
% %                         end

                  end
            end
%   for nn=1:N 
%     for mm=nn+1:N    
%        dy=y(mm)-y(nn);            
%        if(abs(dy)<D) 
%           dx=x(mm)-x(nn);                      
%           dnm=dx.^2+dy.^2;          
%           if(dnm<D^2)                
%             dnm=sqrt(dnm);           
%             F=-(D/dnm);     
% 
%             Fx(nn)= F.*dx;     
%             Fx(mm)= - F.*dx;     
%             Fy(nn)= F.*dy;    
%             Fy(mm)=- F.*dy; 
% 
% %             Fx(nn)=Fx(nn) + F.*dx;     
% %             Fx(mm)=Fx(mm) - F.*dx;     
% %             Fy(nn)=Fy(nn) + F.*dy;    
% %             Fy(mm)=Fy(mm) - F.*dy;
% 
%           end
%         end
%      end
%  end
%             end   
        end  
% 
%  for nn=1:N 
%     for mm=nn+1:N    
%        dy=y(mm)-y(nn);            
%        if(abs(dy)<D) 
%           dx=x(mm)-x(nn);                      
%           dnm=dx.^2+dy.^2;          
%           if(dnm<D^2)                
%             dnm=sqrt(dnm);           
%             F=-(D/dnm);     
% 
%             Fx(nn)= F.*dx;     
%             Fx(mm)= - F.*dx;     
%             Fy(nn)= F.*dy;    
%             Fy(mm)=- F.*dy; 
% 
% %             Fx(nn)=Fx(nn) + F.*dx;     
% %             Fx(mm)=Fx(mm) - F.*dx;     
% %             Fy(nn)=Fy(nn) + F.*dy;    
% %             Fy(mm)=Fy(mm) - F.*dy;
% 
%           end
%         end
%      end
%  end

aax=Fx./M;                
ay=Fy./M;

vx=vx+(ax_old+aax)*dt/2;   
vy=vy+(ay_old+ay)*dt/2;

vxs(i,:)=vx;                 
vys(i,:)=vy;

aax_old=aax;               
ay_old=ay;

plot(t(i),p1,'co','MarkerSize', p1size, 'MarkerFaceColor','c');%%%%NFT%%%
hold on
scatter (xs(i,:),ys(i,:), 'filled')
plot(t(i),p1,'co','MarkerSize', p1size, 'MarkerFaceColor','c'); % plots location of blue particle
ax.NextPlot = 'replacechildren'; % plot replaces location of ribbon, without resetting axes properties
% captures plot as a frame in movie vector
mov(i) = getframe(anim);
    
end
close

% create "object" name "movie" to write our video file to 
movie = VideoWriter(moviename, 'MPEG-4');
movie.FrameRate = fps;
% Write frames collected in mov to "movie" 
open(movie);
writeVideo(movie,mov);
close(movie); 
