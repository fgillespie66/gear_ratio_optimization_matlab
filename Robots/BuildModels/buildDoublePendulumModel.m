function model = buildDoublePendulumModel(varargin)
% this function generates a double pendulum model the user has two options:
        % 1) pass a structure of model parameters with fields 
                % l1_mass
                % l1_l
                % l1_r
                % l1_cm
                % l1_I
                % l2_mass
                % l2_l
                % l2_r
                % l2_cm
                % l2_I
        % 2) no argument : a set of default model parameters are generated
        % using the function get_default_params() defined below

% Model Parameters
if isempty(varargin)
    params = get_default_params(); % should be passed argument to this function
else
    params = varargin{1};
end

model.params = params;

% Actuation
model.actuated_dofs = [1 2];
model.num_act = length(model.actuated_dofs);


% Initialize model struct:
model.NB = 2;
model.gravity = [0 0 -9.81];                   % gravity
model.parent  = zeros(1,model.NB);             % parent body indices
model.jtype   = repmat({'  '},model.NB,1);     % joint types
model.Xtree   = repmat({eye(6)},model.NB,1);   % coordinate transforms
model.I       = repmat({zeros(6)},model.NB,1); % spatial inertias
nb = 0; % current body index

% Pendulum 1
nb = nb + 1;
model.parent(nb) = nb - 1;
model.jtype{nb}  = 'Ry';
model.Xtree{nb}  = eye(6);
model.I{nb}      = mcI(params.l1_mass,params.l1_cm,...
    params.l1_I);
cyl1 = [0 0 0;...
    0 0 params.l1_l];
model.appearance.body{nb} = ...
    {'colour',[245 183 76]./255,...
    'cyl',cyl1,params.l1_r};

% Pendulum 2
nb = nb + 1;
model.parent(nb) = nb - 1;
model.jtype{nb}  = 'Ry';
model.Xtree{nb}  = plux(eye(3), [0 0 params.l1_l]');
model.I{nb}      = mcI(params.l2_mass,params.l2_cm,...
    params.l2_I);
cyl2 = [0 0 0;0 0 params.l2_l];
cylJoint = [0 -1.1*params.l1_r 0;0 1.1*params.l1_r 0];
model.appearance.body{nb} = ...
    {'colour',[89 205 247]./255,...
    'cyl', cyl2, params.l2_r,...
    'cyl',cylJoint,params.l2_r};

% Rail / Background
cylJoint = [0 -1.1*params.l1_r 0;0 1.1*params.l1_r 0];
model.appearance.base = ...
    {'colour',[1 1 1],...
    'cyl', cylJoint, params.l1_r,...
    'tiles', [-10 10; 1 1; -1.0 0.5], 0.4};

% Camera
model.camera.direction = [0 -1 0.2];
model.camera.zoom = 0.4;

end

function params = get_default_params()
        % Mass Parameters
        params.l1_mass = 1.5;
        params.l1_l = 0.15;
        params.l1_cm = 0.5*[0 0 params.l1_l]';
        params.l1_r = 0.008;
        params.l1_I =  rodInertia(params.l1_mass,params.l1_r,params.l1_l);
        
        params.l2_mass = 1.5;
        params.l2_l = 0.15;
        params.l2_cm = 0.5*[0 0 params.l2_l]';
        params.l2_r = 0.008;
        params.l2_I =  rodInertia(params.l2_mass,params.l2_r,params.l2_l);
end