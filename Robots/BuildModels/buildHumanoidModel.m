function model = buildHumanoidModel(params)
% Builds the humanoid robot described by the design
% parameters of the input "params"
%
% ************     3D Humanoid Model:           ************
% Coordinate system: x forward, y left, z up
% Link Indices:
% 0   - fixed origin
% 1   - Px (translate fore/aft)   (massless)
% 2   - Py (translate lateral)    (massless)
% 3   - Pz (translate vertical)   (massless)
% 4   - Rz (yaw)                  (massless)
% 5   - Ry (pitch)                (massless)
% 6   - Rx (roll)                 (base)
% 7,8,9,10,11    - Right leg (hipRz,hipRy,hipRx,knee,ankle)
% 12,13,14,15,16 - Left leg  (hipRz,hipRy,hipRx,knee,ankle)
% 17,18,19,20    - Right arm (shoulderRy,shoulderRx,shoulderRz,elbow)
% 21,22,23,24    - Left arm  (shoulderRy,shoulderRx,shoulderRz,elbow)
%
% Ground contact points:
% 1 - Right toe
% 2 - Left toe
% 3 - Right heel
% 4 - Left heel

%% Various mass parameters and locations
I_torso      = mcI(params.torso_mass, params.torso_COM, reshape(params.torsoRotationalInertia,3,3));
I_shoulderRy = mcI(params.shoulderRyMass,params.shoulderRyCOM,reshape(params.shoulderRyRotationalInertia,3,3));
I_shoulderRx = mcI(params.shoulderRxMass,params.shoulderRxCOM,reshape(params.shoulderRxRotationalInertia,3,3));
I_shoulderRz = mcI(params.shoulderRzMass,params.shoulderRzCOM,reshape(params.shoulderRzRotationalInertia,3,3));
I_elbow      = mcI(params.elbowMass,params.elbowCOM,reshape(params.elbowRotationalInertia,3,3));
I_hipRz      = mcI(params.hipRzMass,params.hipRzCOM,reshape(params.hipRzRotationalInertia,3,3));
I_hipRx      = mcI(params.hipRxMass,params.hipRxCOM,reshape(params.hipRxRotationalInertia,3,3));
I_hipRy      = mcI(params.hipRyMass,params.hipRyCOM,reshape(params.hipRyRotationalInertia,3,3));
I_knee       = mcI(params.kneeMass,params.kneeCOM,reshape(params.kneeRotationalInertia,3,3));
I_ankle      = mcI(params.ankleMass,params.ankleCOM,reshape(params.ankleRotationalInertia,3,3));

%% Initialize model struct:
NLEGS = 2;
NARMS = 2;
N_GC  = 2*NLEGS;
model.NB = 24;                                     % number of bodies
model.NLEGS = NLEGS;                               % number of legs
model.NARMS = NARMS;                               % number of arms
model.N_GC  = N_GC;                                % number of ground contact points
model.gravity   = [0 0 -9.81];                     % gravity
model.parent    = zeros(1,model.NB);               % parent body indices
model.jtype     = repmat({'  '},model.NB,1);       % joint types
model.Xtree     = repmat({eye(6)},model.NB,1);     % coordinate transforms
model.I         = repmat({zeros(6)},model.NB,1);   % spatial inertias
model.gc_X      = repmat({zeros(6)},model.N_GC,1); % coordinate Xforms to gc
model.gc_loc    = repmat({zeros(3)},model.N_GC,1); % position of gc relative to parent frame
model.gc_parent = zeros(1,model.N_GC);

% Define a "home" state of the robot
q_arm = [0 0 0 -0.1];
q_leg = deg2rad([0 0 -30 60 -30]);
model.q_home = [zeros(1,6) repmat(q_leg,1,NLEGS) repmat(q_arm,1,NARMS)];

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
model.jtype{nb}  = 'Rz';
model.Xtree{nb}  = eye(6);
model.I{nb}      = zeros(6);
% Base pitch (no mass)
nb = nb + 1;
model.parent(nb) = nb - 1;
model.jtype{nb}  = 'Ry';
model.Xtree{nb}  = eye(6);
model.I{nb}      = zeros(6);
% Base roll (body mass)
nb = nb + 1;
model.parent(nb) = nb - 1;
model.jtype{nb}  = 'Rx';
model.Xtree{nb}  = eye(6);
model.I{nb}      = I_torso;

% Loop through legs
side_sign = [1 1;1 -1;1 1];
nb_base = nb;
for leg = 1:NLEGS
    % Hip Rz
    nb = nb + 1;
    model.parent(nb) = nb_base;
    model.jtype{nb}  = 'Rz';
    XRot = ry(params.hipRzPitch);
    model.Xtree{nb}  = plux(XRot,side_sign(:,leg)'.*params.hipRzLocation);
    if leg == 2
        model.I{nb}      = I_hipRz;
    else
        model.I{nb}      = flipAlongAxis(I_hipRz,'Y');
    end
    
    % Hip Rx
    nb = nb + 1;
    model.parent(nb) = nb - 1;
    model.jtype{nb}  = 'Rx';
    XRot = ry(params.hipRxPitch);
    model.Xtree{nb}  = plux(XRot,side_sign(:,leg)'.*params.hipRxLocation);
    if leg == 2
        model.I{nb}      = I_hipRx;
    else
        model.I{nb}      = flipAlongAxis(I_hipRx,'Y');
    end
    
    % Hip Ry
    nb = nb + 1;
    model.parent(nb) = nb - 1;
    model.jtype{nb}  = 'Ry';
    XRot = ry(params.hipRyPitch);
    model.Xtree{nb}  = plux(XRot,side_sign(:,leg)'.*params.hipRyLocation);
    if leg == 2
        model.I{nb}      = I_hipRy;
    else
        model.I{nb}      = flipAlongAxis(I_hipRy,'Y');
    end
    
    % Knee
    nb = nb + 1;
    model.parent(nb) = nb - 1;
    model.jtype{nb}  = 'Ry';
    model.Xtree{nb}  = plux(eye(3),side_sign(:,leg)'.*params.kneeLocation);
    if leg == 2
        model.I{nb}      = I_knee;
    else
        model.I{nb}      = flipAlongAxis(I_knee,'Y');
    end
    
    % Ankle
    nb = nb + 1;
    model.parent(nb) = nb - 1;
    model.jtype{nb}  = 'Ry';
    model.Xtree{nb}  = plux(eye(3),side_sign(:,leg)'.*params.ankleLocation);
    if leg == 2
        model.I{nb}      = I_ankle;
    else
        model.I{nb}      = flipAlongAxis(I_ankle,'Y');
    end
       
end

% Loop through arms
for arm = 1:NARMS
    % Shoulder Ry
    nb = nb + 1;
    model.parent(nb) = nb_base;
    model.jtype{nb}  = 'Ry';
    model.Xtree{nb}  = plux(eye(3),side_sign(:,arm)'.*params.shoulderRyLocation);
    if arm == 2
        model.I{nb}      = I_shoulderRy;
    else
        model.I{nb}      = flipAlongAxis(I_shoulderRy,'Y');
    end
    
    % Shoulder Rx
    nb = nb + 1;
    model.parent(nb) = nb - 1;
    model.jtype{nb}  = 'Rx';
    model.Xtree{nb}  = plux(eye(3),side_sign(:,arm)'.*params.shoulderRxLocation);
    if arm == 2
        model.I{nb}      = I_shoulderRx;
    else
        model.I{nb}      = flipAlongAxis(I_shoulderRx,'Y');
    end
    
    % Shoulder Rz
    nb = nb + 1;
    model.parent(nb) = nb - 1;
    model.jtype{nb}  = 'Rz';
    model.Xtree{nb}  = plux(eye(3),side_sign(:,arm)'.*params.shoulderRzLocation);
    if arm == 2
        model.I{nb}      = I_shoulderRz;
    else
        model.I{nb}      = flipAlongAxis(I_shoulderRz,'Y');
    end

    
    % Elbow
    nb = nb + 1;
    model.parent(nb) = nb - 1;
    model.jtype{nb}  = 'Ry';
    model.Xtree{nb}  = plux(eye(3),side_sign(:,arm)'.*params.elbowLocation);
    if arm == 2
        model.I{nb}      = I_elbow;
    else
        model.I{nb}      = flipAlongAxis(I_elbow,'Y');
    end
    
end

%% index selection
%floating base
model.indx.xyz = [1 2 3]; % torso xyz position
% model.indx.x = [1 4];
% model.indx.y = [2 5];
% model.indx.z = [3 6];
model.indx.rpy = [4 5 6]; % torso rpy

%legs
model.indx.hipz = [7 12];
model.indx.hipx = [8 13];
model.indx.hipy = [9 14];
model.indx.knee = [10 15];
model.indx.ankle = [11 16];
%arms
model.indx.shouldery = [17 21];
model.indx.shoulderx = [18 22];
model.indx.shoulderz = [19 23];
model.indx.elbow = [20 24];

model.indx.leg = sort([model.indx.hipx model.indx.hipy model.indx.hipz model.indx.knee model.indx.ankle]);
model.indx.arm = sort([model.indx.shouldery model.indx.shoulderx model.indx.shoulderz model.indx.elbow]);

model.n_leg_joints = 5;
model.n_arm_joints = 4;

%% Ground Contact Points
n_gc = 0;
% Right toe
n_gc = n_gc + 1;
model.gc_parent(n_gc) = model.indx.ankle(1);
model.gc_X{n_gc}      = plux(eye(3),[params.footToeLength 0 -params.footHeight]');
model.gc_loc{n_gc}    = [params.footToeLength 0 -params.footHeight];
% Left toe
n_gc = n_gc + 1;
model.gc_parent(n_gc) = model.indx.ankle(2);
model.gc_X{n_gc}      = plux(eye(3),[params.footToeLength 0 -params.footHeight]');
model.gc_loc{n_gc}    = [params.footToeLength 0 -params.footHeight];
% Right heel
n_gc = n_gc + 1;
model.gc_parent(n_gc) = model.indx.ankle(1);
model.gc_X{n_gc}      = plux(eye(3),[-params.footHeelLength 0 -params.footHeight]');
model.gc_loc{n_gc}    = [-params.footHeelLength 0 -params.footHeight];
% Left heel
n_gc = n_gc + 1;
model.gc_parent(n_gc) = model.indx.ankle(2);
model.gc_X{n_gc}      = plux(eye(3),[-params.footHeelLength 0 -params.footHeight]');
model.gc_loc{n_gc}    = [-params.footHeelLength 0 -params.footHeight];

% Right Hand
n_gc = n_gc + 1;
model.gc_parent(n_gc) = model.indx.elbow(1);
model.gc_X{n_gc}      = plux(eye(3),[0 0 -params.lowerArmLength]');
model.gc_loc{n_gc}    = [0 0 -params.lowerArmLength];

% Left Hand
n_gc = n_gc + 1;
model.gc_parent(n_gc) = model.indx.elbow(2);
model.gc_X{n_gc}      = plux(eye(3),[0 0 -params.lowerArmLength]');
model.gc_loc{n_gc}    = [0 0 -params.lowerArmLength];
        
model.gc_indx.toe  = [1 2];
model.gc_indx.heel = [3 4];
model.gc_indx.hand = [5 6];

%% Actuators (leg + arm)
model.gr = [params.hipGearRatio;params.hipGearRatio;...
    params.hipGearRatio;params.kneeGearRatio;params.ankleGearRatio;...
    params.shoulderGearRatio;params.shoulderGearRatio;...
    params.shoulderGearRatio;params.elbowGearRatio];
model.kt = [params.smallMotorKT;params.smallMotorKT;...
    params.largeMotorKT;params.largeMotorKT;params.smallMotorKT;...
    params.smallMotorKT;params.smallMotorKT;...
    params.smallMotorKT;params.smallMotorKT];
model.Rm = [params.smallMotorR;params.smallMotorR;...
    params.largeMotorR;params.largeMotorR;params.smallMotorR;...
    params.smallMotorR;params.smallMotorR;...
    params.smallMotorR;params.smallMotorR];
tauMax = model.gr.*[params.smallMotorTauMax;params.smallMotorTauMax;...
    params.largeMotorTauMax;params.largeMotorTauMax;params.smallMotorTauMax;...
    params.smallMotorTauMax;params.smallMotorTauMax;...
    params.smallMotorTauMax;params.smallMotorTauMax];
%model.tauMax = [tauMax;tauMax];
model.tauMax = [tauMax(1:5);tauMax(1:5);tauMax(6:9);tauMax(6:9)];
model.batteryV = params.batteryV;

model.ts_int    = repmat([40 40 92 200 70]',2,1);
model.ts_slope  = repmat([0.72 0.72 2 9 2]',2,1);
model.motor_mod = repmat([2 2 3 3 2]',2,1);
model.belt      = repmat([1 1 1 2 1.555]',2,1);

%% Joint Limits
% Joint Limits
hipY_lim  = [-90 90];
knee_lim  = [0 170];
ankle_lim = [-150 150];
right_leg_lim = [-150 75;-150 75;hipY_lim;knee_lim;ankle_lim];
left_leg_lim  = [-75 150;-75 150;hipY_lim;knee_lim;ankle_lim];

shoulderY_lim = [-120 120];
shoulderZ_lim = [-180 180];
elbow_lim     = [-170 0];
right_arm_lim = [shoulderY_lim;-150 20;shoulderZ_lim;elbow_lim]; 
left_arm_lim  = [shoulderY_lim;-20 150;shoulderZ_lim;elbow_lim]; 

model.joint_lim = [right_leg_lim;left_leg_lim;right_arm_lim;left_arm_lim];
model.joint_lim = deg2rad(model.joint_lim);

% Joint velocity limits
leg_joint_vel_limit = [params.smallMotorVelMax;params.smallMotorVelMax;...
    params.largeMotorVelMax;params.largeMotorVelMax;params.smallMotorVelMax]./model.gr(1:5);
arm_joint_vel_limit = [params.smallMotorVelMax;params.smallMotorVelMax;...
    params.smallMotorVelMax;params.smallMotorVelMax]./model.gr(6:9);
model.joint_vel_lim = [repmat(leg_joint_vel_limit,2,1);repmat(arm_joint_vel_limit,2,1)];

%% Model Mass
model.fb_type = 'eul';
H = get_mass_matrix(model, zeros(model.NB,1));
model.mass = H(6,6);

disp(['Created Humanoid model with ' num2str(nb) ' coordinates!']);
