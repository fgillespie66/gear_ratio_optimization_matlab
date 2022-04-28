function model = buildCartPoleModel(params)
% Builds a simple cart pole model (no rotor inertia)
%
% TODO: make params input optional - if missing, then have some default
% parameters

%% Initialize Model Struct
model.NB = 2;
model.gravity = [0 0 -9.81];
model.parent  = zeros(1,model.NB);
model.jtype      = repmat({' '},model.NB,1);       % joint types
model.Xtree      = repmat({eye(6)},model.NB,1);     % coordinate Xforms
model.I          = repmat({zeros(6)},model.NB,1);   % spatial inertias

%% Build the model
% Cart
nb = 1;
model.parent(nb) = nb - 1;
model.jtype{nb}  = 'Px';
model.Xtree{nb}  = eye(6);
a_cart = zeros(10,1);
a_cart(1) = params.cart.mass;
model.I{nb} = inertiaVecToMat(a_cart);
model.appearance.body{nb} = ...
    {'colour',[0.8 0.2 0.2],...
    'box', [-0.5*params.cart.width -0.5*params.cart.depth -0.5*params.cart.height;...
    0.5*params.cart.width 0.5*params.cart.depth 0.5*params.cart.height]};

% Pole
nb = 2;
model.parent(nb) = nb - 1;
model.jtype{nb}  = 'Ry';
model.Xtree{nb}  = rotz(pi);
a_pole = zeros(10,1);
a_pole(1) = params.pole.mass;
a_pole(4) = 0.5*params.pole.length;
a_pole(6) = params.pole.Iyy;
model.I{nb} = inertiaVecToMat(a_pole);
model.appearance.body{nb} = ...
    {'colour',[0.2 0.4 0.78],...
    'cyl', [0 0 0;0 0 params.pole.length], params.pole.rad};

% Base
model.appearance.base = ...
    {'colour',[0.91 0.91 0.91],...
    'cyl', [-50 0 0;50 0 0], params.base.rad,...
    'tiles',[-50 50;0.5 0.5;-20 20],0.4};

%% Camera
model.camera.body = 1;
model.camera.direction = [0.05 -0.6 0.2];
model.camera.zoom = 0.65;

disp(['Created Cart Pole model with ' num2str(nb) ' coordinates!']);

