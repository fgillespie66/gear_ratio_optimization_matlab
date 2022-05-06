function animate_simple(t,z,kinematics,speed, gear_ratio, body_width)

    axis([-.4 .5 -.1 1.2])
    h_ground = plot([-1 1],[0 0],'k-','LineWidth',5);
    hold on
    title("Gear Ratio: " + gear_ratio)
    h_leg_l    = plot([0],[0],'-o',...
                'LineWidth',3,..., 
                'MarkerEdgeColor','r',...
                'MarkerFaceColor','r',...
                'MarkerSize',6); 

    %h_leg_r    = plot([0],[0],'-o',...
    %            'LineWidth',3,..., 
    %            'MarkerEdgeColor','r',...
    %            'MarkerFaceColor','r',...
    %            'MarkerSize',6); 

    %body    = plot([0],[0],'-o',...
    %            'LineWidth',3,..., 
    %            'MarkerEdgeColor','r',...
    %            'MarkerFaceColor','r',...
    %            'MarkerSize',6); 

    tic                                             % start counter
    while toc < t(end)/speed                        % while there's real time left
        tsim = toc*speed;                           % determine the simulation time to draw
        zint = interp1(t',z',tsim', 'linear')';     % interpolate to get coordinates at that time
        z1 = draw_lines(zint,kinematics,h_leg_l, 0);
        %z1 = draw_lines(zint,kinematics,h_leg_l, body_width);
        %z2 = draw_lines(zint,kinematics,h_leg_r, -body_width);
        %body.XData = [-body_width/2, body_width/2];
        %body.YData = [z1, z2];
        drw_now
    end
    z1end = draw_lines(zint,kinematics,h_leg_l, 0);
    %z1end = draw_lines(z(:,end),kinematics,h_leg_l, body_width);
    %z2end = draw_lines(z(:,end),kinematics,h_leg_r, -body_width);
    %body.XData = [-body_width/2, body_width/2];
    %body.YData = [z1end, z2end];
    drw_now;
end

function zh = draw_lines(z,kinematics,h_leg_l, body_width)
    keypoints = full( kinematics.keypoints(z) );
    h_leg_l.XData = keypoints(1,:) + body_width/2;
    h_leg_l.YData = keypoints(2,:);
    zh = h_leg_l.YData(1);
end

function drw_now()
    drawnow;
    axis equal % sets the X:Y aspect ratio 1:1; otherwise things will look stretched
    axis([-.4 .5 -.1 1.2])
end
