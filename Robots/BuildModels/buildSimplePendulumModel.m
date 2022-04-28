function model = buildSimplePendulumModel(varargin)
% this function generates a simple pendulum model the user has two options:
        % 1) pass a structure of model parameters with fields 
                % l1_mass
                % l1_l
                % l1_r
                % l1_cm
                % l1_I

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
model.actuated_dofs = [1];
model.num_act = length(model.actuated_dofs);
model.u_lim = params.u_lim;


% Initialize model struct:
model.NB = 1;
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
        params.l1_mass = 2.0;
        params.l1_l = 0.4;
        params.l1_cm = 0.5*[0 0 params.l1_l]';
        params.l1_r = 0.008;
        params.l1_I =  rodInertia(params.l1_mass,params.l1_r,params.l1_l);
        params.u_lim = 2.5;
end