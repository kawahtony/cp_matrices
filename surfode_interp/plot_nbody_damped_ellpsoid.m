clear ; clc ; close all ;

% azure = [0, 128, 255]/255 ;
jade = [0, 168, 107]/255 ;


%% Load file

numpt = 5 ;
filename = ['nbody_ellipsoid_damped_num=', num2str(numpt), '.mat'] ;
load(filename) ;


%% Video setup

videoName = ['nbody_ellipsoid_damped_num=', num2str(numpt), '.mp4'] ;
vidObj = VideoWriter(videoName, 'MPEG-4') ;
vidObj.FrameRate = 60 ;
open(vidObj) ;


%% Initialize plot


figure(1) ; clf ; hold on ;
ellipsoid(0, 0, 0, 1, 1, 1.5, 50, FaceColor = jade, FaceAlpha = 0.5, EdgeColor = 'none')
axis equal
view(3)

% distinct color per particle, shared by its marker and its trail
colors = lines(numpt) ;

x = x_stored(:,:,1) ;
m1 = scatter3(x(1,:), x(2,:), x(3,:), 36, colors, 'filled', 'MarkerEdgeColor', 'k') ;
titleHandle = title('t = 0', 'FontSize', 16) ;

% sliding window of the most recent trail points (fixed length, oldest
% points drop off as new ones are recorded)
trailLen = 50 ;
trail_x = nan(numpt, trailLen) ;
trail_y = nan(numpt, trailLen) ;
trail_z = nan(numpt, trailLen) ;

trail_x(:,end) = x(1,:)' ;
trail_y(:,end) = x(2,:)' ;
trail_z(:,end) = x(3,:)' ;

% one trail line per particle, colored to match its marker
trailHandles = gobjects(numpt,1) ;
for k = 1:numpt
    trailHandles(k) = plot3(trail_x(k,:), trail_y(k,:), trail_z(k,:), ...
        '-', 'Color', colors(k,:), 'LineWidth', 1) ;
end

% capture the first frame
drawnow() ;
writeVideo(vidObj, getframe(gcf)) ;

%% Looping over time


for i = 1:size(x_stored, 3)

    x = x_stored(:,:,i) ;
    t = t_stored(i) ;

    % update marker positions in place (no delete/recreate)
    set(m1, 'XData', x(1,:), 'YData', x(2,:), 'ZData', x(3,:)) ;

    % slide the trail window forward and redraw for each particle
    trail_x = [trail_x(:,2:end), x(1,:)'] ;
    trail_y = [trail_y(:,2:end), x(2,:)'] ;
    trail_z = [trail_z(:,2:end), x(3,:)'] ;
    for k = 1:numpt
        set(trailHandles(k), ...
            'XData', trail_x(k,:), ...
            'YData', trail_y(k,:), ...
            'ZData', trail_z(k,:)) ;
    end

    drawnow();
    set(titleHandle, 'String', ['t = ' num2str(round(t,1))]) ;

    % grab the current figure and append it as a video frame
    writeVideo(vidObj, getframe(gcf)) ;

end

close(vidObj) ;


















% clear ; clc ; close all ;
% 
% % azure = [0, 128, 255]/255 ;
% jade = [0, 168, 107]/255 ;
% 
% 
% %% Load file
% 
% numpt = 5 ;
% filename = ['nbody_ellipsoid_damped_num=', num2str(numpt), '.mat'] ;
% load(filename) ;
% 
% 
% %% Initialize plot
% 
% 
% figure(1) ; clf ; hold on ;
% ellipsoid(0, 0, 0, 1, 1, 1.5, 50, FaceColor = jade, FaceAlpha = 0.5, EdgeColor = 'none')
% axis equal
% view(3)
% 
% % distinct color per particle, shared by its marker and its trail
% colors = lines(numpt) ;
% 
% x = x_stored(:,:,1) ;
% m1 = scatter3(x(1,:), x(2,:), x(3,:), 36, colors, 'filled', 'MarkerEdgeColor', 'k') ;
% titleHandle = title('t = 0', 'FontSize', 16) ;
% 
% % sliding window of the most recent trail points (fixed length, oldest
% % points drop off as new ones are recorded)
% trailLen = 50 ;
% trail_x = nan(numpt, trailLen) ;
% trail_y = nan(numpt, trailLen) ;
% trail_z = nan(numpt, trailLen) ;
% 
% trail_x(:,end) = x(1,:)' ;
% trail_y(:,end) = x(2,:)' ;
% trail_z(:,end) = x(3,:)' ;
% 
% % one trail line per particle, colored to match its marker
% trailHandles = gobjects(numpt,1) ;
% for k = 1:numpt
%     trailHandles(k) = plot3(trail_x(k,:), trail_y(k,:), trail_z(k,:), ...
%         '-', 'Color', colors(k,:), 'LineWidth', 1) ;
% end
% 
% %% Looping over time
% 
% 
% for i = 1:size(x_stored, 3)
% 
%     x = x_stored(:,:,i) ;
%     t = t_stored(i) ;
% 
%     % update marker positions in place (no delete/recreate)
%     set(m1, 'XData', x(1,:), 'YData', x(2,:), 'ZData', x(3,:)) ;
% 
%     % slide the trail window forward and redraw for each particle
%     trail_x = [trail_x(:,2:end), x(1,:)'] ;
%     trail_y = [trail_y(:,2:end), x(2,:)'] ;
%     trail_z = [trail_z(:,2:end), x(3,:)'] ;
%     for k = 1:numpt
%         set(trailHandles(k), ...
%             'XData', trail_x(k,:), ...
%             'YData', trail_y(k,:), ...
%             'ZData', trail_z(k,:)) ;
%     end
% 
%     drawnow();
%     set(titleHandle, 'String', ['t = ' num2str(round(t,1))]) ;
% 
% end