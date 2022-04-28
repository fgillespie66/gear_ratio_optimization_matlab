function [gh,gg_x,gg_y]  = getGroundInfo(terrain, pf)
% Outputs
% gh = ground height
% gg = ground gradients in the x and y directions

switch terrain.type
    case 'sigmoid3D'
        a = terrain.z_height; 
        b = terrain.x_shift;
        c = terrain.steepness;
        d = terrain.y_angle;
        
        gh = a/(1+exp(-c*(pf(1)-b-d*pf(2))));
        gg_x = ( a*c*exp(-c*(pf(1)-d*pf(2)-b)) ) / ( (1+exp(-c*(pf(1)-b-d*pf(2))))^2 );
        gg_y = ( a*c*(-d)*exp(-c*(pf(1)-d*pf(2)-b)) ) / ( (1+exp(-c*(pf(1)-b-d*pf(2))))^2 );
        
    case 'sigmoid3D_2ndOrderY'
        a = terrain.z_height; 
        b = terrain.x_shift;
        c = terrain.steepness;
        d = terrain.y_angle;
        f = terrain.y_curvature;
        
        gh = a/(1+exp(-c*(pf(1)-b-d*pf(2)-f*pf(2)^2)));
        gg_x = ( a*c*exp(-c*(pf(1)-f*pf(2)^2-d*pf(2)-b)) ) / ( (1+exp(-c*(pf(1)-b-d*pf(2)-f*pf(2)^2)))^2 );
        gg_y = ( a*c*(-2*f*pf(2)-d)*exp(-c*(pf(1)-f*pf(2)^2-d*pf(2)-b)) ) / ( (1+exp(-c*(pf(1)-b-d*pf(2)-f*pf(2)^2)))^2 );
        
    case 'sigmoid2D'
        a = terrain.z_height;
        b = terrain.x_shift;
        c = terrain.steepness;
        
        gh = a/(1+exp(-c*(pf(1)-b)));
        gg_x = ( a*c*exp(-c*(pf(1)-b)) ) / ( (1+exp(-c*(pf(1)-b)))^2 );
        
    case 'obstacle'
        a = terrain.obstacle.height;
        ep = terrain.obstacle.width;
        c = terrain.obstacle.location_x;
        
        gh = a*exp(-(ep.*(pf(1)-c)).^2);
        gg_x = -2*a*ep*(pf(1)-c)*exp(-(ep.*(pf(1)-c)).^2);
        
    case 'obstacle2'
        a = terrain.obstacle.height;
        ep = terrain.obstacle.width;
        c = terrain.obstacle.location_x;
        
        gh = a*exp(-(ep.*(pf(1)-c)).^2);
        gg_x = -2*a*ep*(pf(1)-c)*exp(-(ep.*(pf(1)-c)).^2);
        
    otherwise
        error(['The terrain type ',terrain.type,' is not permissible.'])
        
end


end