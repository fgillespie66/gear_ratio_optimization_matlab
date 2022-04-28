function model = buildQuadrupedModel2D(params)
% Builds the quadruped robot described by the design
% parameters of the input "params"
%
% ************     2D Quadruped Model:           ************
% Coordinate system: x forward, y up
% Link Indices:
% 0   - fixed origin
% 1   - px (translate right/left) (massless)
% 2   - py (translate vertical)   (massless)
% 3   - r  (body pitch)           (base)
% 4,5 - Front leg                 (hip,knee)
% 6,7 - Rear leg                  (hip,knee)
% 
% Ground contact points:
% 1 - front foot
% 2 - rear foot

%% Initialize model struct:
model.NB         = 7;                               % number of bodies
model.NLEGS      = 2;                               % number of legs
model.N_GC       = model.NLEGS;                     % number of ground contacts
model.NARMS      = 0;                               % number of arms
model.gravity    = [0 -9.81];                       % gravity
model.parent     = zeros(1,model.NB);               % parent body indices
model.jtype      = repmat({'  '},model.NB,1);       % joint types
model.Xtree      = repmat({eye(3)},model.NB,1);     % coordinate Xforms
model.I          = repmat({zeros(3)},model.NB,1);   % spatial inertias
model.gc_X       = repmat({zeros(3)},model.N_GC,1); % coordinate Xforms to gc
model.gc_loc     = repmat({zeros(1,2)},model.N_GC,1); % position of gc relative to parent frame
model.gc_parent  = zeros(1,model.N_GC);

% Define a "home" position of the robot
q_leg = [-1.15 2.25];
model.q_home  = [zeros(1,3) repmat(q_leg,1,model.NLEGS)];

%% Build the Model
nb = 0; % current body index
% Base x translation (no mass)
nb = nb + 1;
model.parent(nb) = nb - 1;
model.jtype{nb}  = 'px';
model.Xtree{nb}  = eye(3);
model.I{nb}      = zeros(3);
% Base y translation (no mass)
nb = nb + 1;
model.parent(nb) = nb - 1;
model.jtype{nb}  = 'py';
model.Xtree{nb}  = eye(3); 
model.I{nb}      = zeros(3);
% Base pitch
nb = nb + 1;
model.parent(nb) = nb - 1;
model.jtype{nb} = 'r';
model.Xtree{nb} = eye(3);
model.I{nb}     = params.bodyInertia;

nb_base = nb;
% Front Hip
nb = nb + 1;
model.parent(nb) = nb_base;
model.jtype{nb}  = 'r';
model.Xtree{nb}  = plnr(0, [1 1].*params.abadLocation([1 3]) );
model.I{nb}  = params.hipInertia;

% Front Knee
nb = nb + 1;
model.parent(nb) = nb-1;
model.jtype{nb}  = 'r';
model.Xtree{nb}  = plnr(0, params.kneeLocation([1 3]) );
model.I{nb}  = params.kneeInertia;

% Rear Hip
nb = nb + 1;
model.parent(nb) = nb_base;
model.jtype{nb}  = 'r';
model.Xtree{nb}  = plnr(0, [-1 1].*params.abadLocation([1 3]) );
model.I{nb}  = params.hipInertia;

% Rear Knee
nb = nb + 1;
model.parent(nb) = nb-1;
model.jtype{nb}  = 'r';
model.Xtree{nb}  = plnr(0, params.kneeLocation([1 3]) );
model.I{nb}  = params.kneeInertia;

model.hip_idx  = [4 6];
model.knee_idx = [5 7];

% Front foot
model.gc_X{1}      = plnr(0, params.footLocation([1 3]) );
model.gc_loc{1}    = params.footLocation([1 3]);
model.gc_parent(1) = model.knee_idx(1);

% Rear foot
model.gc_X{2}      = plnr(0, params.footLocation([1 3]) );
model.gc_loc{2}    = params.footLocation([1 3]);
model.gc_parent(2) = model.knee_idx(2);

% Actuators
model.gr = [params.hipGearRatio;params.kneeGearRatio];
model.kt = [params.motorKT;params.motorKT];
model.Rm = [params.motorR;params.motorR];
tauMax = model.gr.*[params.motorTauMax;params.motorTauMax];
model.tauMax = repmat(tauMax,2,1);
model.batteryV = params.batteryV;

% Joint Limits
hip_lim  = [-90 90];
knee_lim = [0 170];
leg_lim = [hip_lim;knee_lim];
model.joint_lim = repmat(leg_lim,2,1);
model.joint_lim = deg2rad(model.joint_lim);

% Joint velocity limits
joint_vel_limit = [params.motorVelMax;params.motorVelMax]./model.gr;
model.joint_vel_lim = repmat(joint_vel_limit,2,1);

% Model Mass
H = HandC_original(model, zeros(model.NB,1), zeros(model.NB,1));
model.mass = H(1,1);

disp(['Created 2D Quadruped model with ' num2str(nb) ' coordinates!']);

