function model = buildQuadrupedModel_Quat(params)
% Builds the quadruped robot described by the design
% parameters of the input "params"
%
% ************     3D Quadruped Model:           ************
% Coordinate system: x forward, y left, z up
% Link Indices:
% 0   - fixed origin
% 1   - Floating Base
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
model.NB         = 13;                              % number of bodies
model.NLEGS      = 4;                               % number of legs
model.N_GC       = model.NLEGS;                     % number of ground contacts
model.NARMS      = 0;                               % number of arms
model.gravity    = [0 0 -9.81];                     % gravity
model.parent     = zeros(1,model.NB);               % parent body indices
model.jtype      = repmat({'  '},1,model.NB);       % joint types
model.Xtree      = repmat({eye(6)},1,model.NB);     % coordinate Xforms
model.I          = repmat({zeros(6)},1,model.NB);   % spatial inertias
model.gc_X       = repmat({zeros(6)},model.N_GC,1); % coordinate Xforms to gc
model.gc_loc     = repmat({zeros(3)},model.N_GC,1); % position of gc relative to parent frame
model.gc_parent  = zeros(1,model.N_GC);

% Define a "home" position of the robot
q_leg = [0 -1.45 2.65];
model.q_home  = [zeros(1,6) repmat(q_leg,1,model.NLEGS)];

%% Build the Model
nb = 0; % current body index
% Floating base
nb = nb + 1;
model.parent(nb) = nb - 1;
model.jtype{nb}  = 'Fb';
model.Xtree{nb}  = eye(6);
model.I{nb}      = params.bodyInertia;

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

model.abad_idx = [2 5 8 11];
model.hip_idx  = [3 6 9 12];
model.knee_idx = [4 7 10 13];

% Loop Through Ground Contacts
for gc = 1:model.N_GC
    model.gc_X{gc}      = plux(eye(3),side_sign(:,leg)'.*params.footLocation);
    model.gc_loc{gc}    = side_sign(:,leg)'.*params.footLocation;
    model.gc_parent(gc) = model.knee_idx(gc);
end

% Actuators
% model.gr = [params.abadGearRatio;params.hipGearRatio;params.kneeGearRatio];
% model.kt = [params.motorKT;params.motorKT;params.motorKT];
% model.Rm = [params.motorR;params.motorR;params.motorR];
% tauMax = model.gr.*[params.motorTauMax;params.motorTauMax;params.motorTauMax];
% model.tauMax = repmat(tauMax,4,1);
% model.batteryV = params.batteryV;

disp(['Created Quadruped model with ' num2str(nb) ' coordinates!']);

end

%% Utility functions
function I = flipAlongAxis(I_in, axis)
h = skew(I_in(1:3,4:6));
Ibar = I_in(1:3,1:3);
m = I_in(6,6);

if strcmp(class(I_in),'casadi.MX') % make sure this function works for optimization parameters
    P = casadi.MX.zeros(4,4);
    I = casadi.MX.eye(6);
else
    P = zeros(4,4);
    I = eye(6);
end

P(1:3,1:3) = 0.5*trace(Ibar)*eye(3) - Ibar;
P(1:3,4) = h;
P(4,1:3) = h';
P(4,4) = m;

X = eye(4);
if (axis == 'X')
    X(1, 1) = -1;
elseif (axis == 'Y')
    X(2, 2) = -1;
elseif (axis == 'Z')
    X(3, 3) = -1;
end
P = X * P * X;

% I = eye(6);
m = P(4,4);
h = P(1:3,4);
E = P(1:3,1:3);
Ibar = trace(E) * eye(3) - E;
I(1:3,1:3) = Ibar;
I(1:3,4:6) = skew(h);
I(4:6,1:3) = skew(h)';
I(4:6,4:6) = m * eye(3);
end
