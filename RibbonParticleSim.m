clear; close all; 
movielength = 30; fps = 60; % length (seconds) and fps of the output viedo
moviename = 'Ribbon';

t =[0:pi/50:2*pi];
y = sin(t);

% particle size 
 p1size = 12;

%setup figure, titles, labels, ...
anim = figure;
title('Ribbon Simulation', 'fontsize', 12); xlabel('x'); ylabel('y');
ax = gca; ax.FontSize = 12;

% plot commands will add lines to an existing figure
ax.NextPlot = 'add';

xlim([-1.2, 10]);
ylim([-1.2, 1.5]);
ax.DataAspectRatio = [1 1 1];

for i = 1:length(t)
    
    plot(t,y)
    hold on
    % Define angle, particle rope at time step
    p1 = sin(t(i));

    plot(t(i),p1,'co','MarkerSize', p1size, 'MarkerFaceColor','c'); % plots location of Ribbon
    ax.NextPlot = 'replacechildren';% plot replaces location of ribbon, without resetting axes properties
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

