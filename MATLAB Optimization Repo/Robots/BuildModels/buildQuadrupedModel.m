function model = buildQuadrupedModel(params)
% Builds the quadruped robot described by the design
% parameters of the input "params"
%
% ************     3D Quadruped Model:           ************
% Coordinate system: x forward, y left, z up
% Link Indices:
% 0   - fixed origin
% 1   - Px (translate fore/aft)   (massless)
% 2   - Py (translate lateral)    (massless)
% 3   - Pz (translate vertical)   (massless)
% 4   - Rz (yaw)                  (massless)
% 5   - Ry (pitch)                (massless)
% 6   - Rx (roll)                 (base)
% 7,8,9    - Front right leg      (ab/ad,hip,knee)
% 10,11,12 - Front left leg       (ab/ad,hip,knee)
% 13,14,15 - Rear right leg       (ab/ad,hip,knee)
% 16,17,18 - Rear left leg        (ab/ad,hip,knee)
% 
% Ground contact points:
% 1 - front right foot
% 2 - front left foot
% 3 - back right foot
% 4 - back left foot

%% Initialize model struct:
model.NB         = 18;                              % number of bodies
model.NLEGS      = 4;                               % number of legs
model.N_GC       = model.NLEGS;                     % number of ground contacts
model.NARMS      = 0;                               % number of arms
model.gravity    = [0 0 -9.81];                     % gravity
model.parent     = zeros(1,model.NB);               % parent body indices
model.jtype      = repmat({'  '},model.NB,1);       % joint types
model.Xtree      = repmat({eye(6)},model.NB,1);     % coordinate Xforms
model.I          = repmat({zeros(6)},model.NB,1);   % spatial inertias
model.gc_X       = repmat({zeros(6)},model.N_GC,1); % coordinate Xforms to gc
model.gc_loc     = repmat({zeros(1,3)},model.N_GC,1); % position of gc relative to parent frame
model.gc_parent  = zeros(1,model.N_GC);

% Define a "home" position of the robot
%q_leg = [0 -1.45 2.65];
q_leg = [0 -1.15 2.25];
model.q_home  = [zeros(1,6) repmat(q_leg,1,model.NLEGS)];

%% Build the Model
nb = 0; % current body index
% Base x translation (no mass)
nb = nb + 1;
model.parent(nb) = nb - 1;
model.jtype{nb}  = 'Px';
model.Xtree{nb}  = eye(6);
model.I{nb}      = zeros(6);
% Base y translation (no mass)
nb = nb + 1;
model.parent(nb) = nb - 1;
model.jtype{nb}  = 'Py';
model.Xtree{nb}  = eye(6); 
model.I{nb}      = zeros(6);
% Base z translation (no mass)
nb = nb + 1;
model.parent(nb) = nb - 1;
model.jtype{nb}  = 'Pz';
model.Xtree{nb}  = eye(6);
model.I{nb}      = zeros(6);
% Base yaw (no mass)
nb = nb + 1;
model.parent(nb) = nb - 1;
model.jtype{nb} = 'Rz';
model.Xtree{nb} = eye(6);
model.I{nb}     = zeros(6);
% Base pitch (no mass)
nb = nb + 1;
model.parent(nb) = nb - 1;
model.jtype{nb} = 'Ry';
model.Xtree{nb} = eye(6);
model.I{nb}     = zeros(6);
% Base roll
nb = nb + 1;
model.parent(nb) = nb - 1;
model.jtype{nb} = 'Rx';
model.Xtree{nb} = eye(6);
model.I{nb}     = params.bodyInertia;

% Loop through legs
nb_base = nb;
side_sign = [1 1 -1 -1;-1 1 -1 1;1 1 1 1];
leg_side = -1;
for leg = 1:model.NLEGS
    % Ab/Ad
    nb = nb + 1;
    model.parent(nb) = nb_base;
    model.jtype{nb}  = 'Rx';
    model.Xtree{nb}  = plux(eye(3),side_sign(:,leg)'.*params.abadLocation);
    if leg_side > 0
        model.I{nb}  = params.abadInertia;
    else
        model.I{nb}  = flipAlongAxis(params.abadInertia,'Y');
    end
    
    % Hip
    nb = nb + 1;
    model.parent(nb) = nb-1;
    model.jtype{nb}  = 'Ry';
    model.Xtree{nb}  = plux(rz(pi),side_sign(:,leg)'.*params.hipLocation);
    if leg_side > 0
        model.I{nb}  = params.hipInertia;
    else
        model.I{nb}  = flipAlongAxis(params.hipInertia,'Y');
    end
    
    % Knee
    nb = nb + 1;
    model.parent(nb) = nb-1;
    model.jtype{nb}  = 'Ry';
    model.Xtree{nb}  = plux(eye(3),side_sign(:,leg)'.*params.kneeLocation);
    if leg_side > 0
        model.I{nb}  = params.kneeInertia;
    else
        model.I{nb}  = flipAlongAxis(params.kneeInertia,'Y');
    end
    
    leg_side = -1 * leg_side;
end

model.abad_idx = [7 10 13 16];
model.hip_idx  = [8 11 14 17];
model.knee_idx = [9 12 15 18];

% Loop Through Ground Contacts
for gc = 1:model.N_GC
    model.gc_X{gc}      = plux(eye(3),side_sign(:,leg)'.*params.footLocation);
    model.gc_loc{gc}    = side_sign(:,leg)'.*params.footLocation;
    model.gc_parent(gc) = model.knee_idx(gc);
end

% Actuators
model.gr = [params.abadGearRatio;params.hipGearRatio;params.kneeGearRatio];
model.kt = [params.motorKT;params.motorKT;params.motorKT];
model.Rm = [params.motorR;params.motorR;params.motorR];
tauMax = model.gr.*[params.motorTauMax;params.motorTauMax;params.motorTauMax];
model.tauMax = repmat(tauMax,4,1);
model.batteryV = params.batteryV;

% Joint Limits
hip_lim  = [-90 90];
knee_lim = [0 170];
right_leg_lim = [-120 80;hip_lim;knee_lim];
left_leg_lim = [-80 120;hip_lim;knee_lim];
model.joint_lim = [right_leg_lim;left_leg_lim;right_leg_lim;left_leg_lim];
model.joint_lim = deg2rad(model.joint_lim);

% Joint velocity limits
joint_vel_limit = [params.motorVelMax;params.motorVelMax;params.motorVelMax]./model.gr;
model.joint_vel_lim = repmat(joint_vel_limit,4,1);

% Model Mass
model.fb_type = 'eul';
H = get_mass_matrix(model, zeros(model.NB,1));
model.mass = H(6,6);

disp(['Created Quadruped model with ' num2str(nb) ' coordinates!']);

