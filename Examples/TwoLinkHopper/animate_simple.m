function animate_simple(t,z,kinematics,speed)

    axis([-.2 .2 -.2 1])
    h_ground = plot([-1 1],[0 0],'k-','LineWidth',5);
    hold on
    h_leg    = plot([0],[0],'-o',...
                'LineWidth',3,...
                'MarkerEdgeColor','r',...
                'MarkerFaceColor','r',...
                'MarkerSize',6); 

    tic                                             % start counter
    while toc < t(end)/speed                        % while there's real time left
        tsim = toc*speed;                           % determine the simulation time to draw
        zint = interp1(t',z',tsim', 'linear')';     % interpolate to get coordinates at that time
        draw_lines(zint,kinematics,h_leg);
    end
    draw_lines(z(:,end),kinematics,h_leg);
end

function draw_lines(z,kinematics,h_leg)
    keypoints = full( kinematics.keypoints(z) );
    h_leg.XData = keypoints(1,:);
    h_leg.YData = keypoints(2,:);
    drawnow
    axis equal % sets the X:Y aspect ratio 1:1; otherwise things will look stretched
    axis([-.2 .2 -.2 1])
end
