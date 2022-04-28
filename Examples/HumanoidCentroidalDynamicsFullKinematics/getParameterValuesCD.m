function pb = getParameterValuesCD(pb)
% shorthand for parameter structure;
p = pb.p;


p.desired_state_cmd.val = [deg2rad(0), 1.5, 0, 0.6]; % [yaw rate, forward velocity (body frame), lateral velocity (body frame), body height] 

% unpack indexes
shouldery_inx = pb.model.indx.shouldery;
shoulderx_inx = pb.model.indx.shoulderx;
shoulderz_inx = pb.model.indx.shoulderz;
elbow_inx = pb.model.indx.elbow;

p.dt.val = pb.dt*ones(pb.N-1,1);

%contact schedule and dt
p.cs.val =  [0     0     0     0     1     1     1     1     0;
             0     0     0     0     1     1     1     1     0;
             1     1     1     1     0     0     0     0     1;
             1     1     1     1     0     0     0     0     1];


% parameters for raibert's heuristic
p.raibert.val = [0.2 0.2]; % raibert gains in xy directions

% COST FUNCTION WEIGHTS
p.Q_step_reg.val = 2*[2000 2000 2000 2000]; %cost function weight.

p.Qq.val = [0 0 1000 2000 2000 2000 0.1.*ones(1,5) ...
    0.1.*ones(1,5) 0.1.*ones(1,8)];

p.Qq.val(pb.model.indx.hipz) =   2000;  % penalize hip Rz  more
p.Qq.val(pb.model.indx.hipx) =  50;  % penalize hip Rx  more

p.Qq.val(shouldery_inx) =   5;  % penalize shoulder Ry  more
p.Qq.val(shoulderx_inx) =   50;  % penalize shoulder Rx  more
p.Qq.val(shoulderz_inx) =   5;  % penalize shoulder Rx  more
p.Qq.val(elbow_inx) =   5;  % penalize shoulder Ry  more


p.Qv.val = [5 5 1000 1000 1000 1000 1.*ones(1,5) ...
    1.*ones(1,5) 0.1.*ones(1,8)];

p.Qv.val(shoulderx_inx) =   0.1;  % penalize shoulder Rx (Rz)  more

p.Qv.val(8) =  50;  % penalize hip Rx  more
p.Qv.val(13) =  50;  % penalize hip Rx  more

p.Qv.val(shouldery_inx(1)) =   0.1;  % penalize shoulder Ry  more
p.Qv.val(shouldery_inx(2)) =   0.1;  % penalize shoulder Ry  more

p.Qr.val = [0 0 0, 0 0 0, 0 0 0];%[rx ry rz,rdx,rdy,rdz,rddx,rddy,rddz]
p.Qh.val = [2000 100 2000, 0 0 0]; % promotes arms swing usage

p.Qc.val = repmat([700 50 700],1,4);
p.Qf.val = repmat([0.01 0.01 0.05],1,4);


% REFERENCE TRAJECTORY
p.pf_ref.val = diag([1 -1 1, 1 1 1, -1 -1 1, -1 1 1])*repmat([0.065 0.07 -0.4],1,pb.model.N_GC)';
%rtoe,ltoe,rheel,lheel

% FORCE BOUNDS
% max normal force
p.fz_max.val = 500;%*ones(pb.model.NLEGS,1); % single contact - [Rtoe;Ltoe;Rheel;Lheel];
p.fz_min.val = 0;%zeros(pb.model.N_GC,1);


% STATE BOUNDS
% Initial state
q_leg = [0 0 -0.4 0.77 -0.37];
q_arm = [0 0 0 -1.0];

q_init = [0 0 0.6565 0 0 0 q_leg q_leg q_arm q_arm]';
qd_init = zeros(pb.model.NB,1);


% Limits
q_init_min = q_init;
q_init_max = q_init;
qd_init_min = qd_init;
qd_init_max = qd_init;

qleg_min = [-pi -pi/2 -deg2rad(120)    0  -pi];
qleg_max = [pi  pi/2  +deg2rad(60) pi pi];
vleg_min = [-35 -35 -40 -50 -50];
vleg_max = -vleg_min;

% CONSERVATIVE ARM JOINT LIMITS to avoid collision
qarm_min = [-1.05*pi,    -pi,        -pi     -deg2rad(150),...
            -1.05*pi   -deg2rad(11)   0      -deg2rad(150)];
        
qarm_max = [1.05*pi     deg2rad(11)   0     +deg2rad(150),...
            1.05*pi,     pi           pi    +deg2rad(150)];

varm_min = [-35 -40 -40 -50]/2;
varm_max = -varm_min;

q_min = [-5000 -5000 0 -deg2rad(1.05*360) -deg2rad(1.1*360) -deg2rad(1.05*360) repmat(qleg_min,1,pb.model.NLEGS) qarm_min]';
q_max = [5000 5000 5 deg2rad(1.05*360) deg2rad(1.1*360) deg2rad(1.05*360) repmat(qleg_max,1,pb.model.NLEGS) qarm_max]';
qd_min = [-99*ones(1,6) repmat(vleg_min,1,pb.model.NLEGS) repmat(varm_min,1,pb.model.NARMS)]';
qd_max = [99*ones(1,6) repmat(vleg_max,1,pb.model.NLEGS) repmat(varm_max,1,pb.model.NARMS)]';

q_final_min = q_min;
q_final_max = q_max;
qd_final_min = qd_min;
qd_final_max = qd_max;

%% task specific bounds

q_min(pb.model.indx.ankle) = [-rad2deg(70) -rad2deg(70)];
q_max(pb.model.indx.ankle) = [pi/2 pi/2];

q_min(pb.model.indx.rpy) = -deg2rad([70 70 8*pi])';
q_max(pb.model.indx.rpy) = deg2rad([70 70 8*pi])';

q_min(pb.model.indx.hipz) = -deg2rad(30);
q_max(pb.model.indx.hipz) = +deg2rad(30);

q_min(pb.model.indx.hipx) = -deg2rad(60);
q_max(pb.model.indx.hipx) = +deg2rad(60);


%% pack problem structure for output

% pack bounds in pb structure
p.q_min.val = q_min;
p.q_max.val = q_max;
p.qd_min.val = qd_min;
p.qd_max.val = qd_max;

% initial state
p.q_init.val  = q_init_min;
p.qd_init.val = qd_init_min;


p.q_final_min.val  = q_final_min;
p.q_final_max.val  = q_final_max;
p.qd_final_min.val = qd_final_min;
p.qd_final_max.val = qd_final_max;

% reference
p.q_des.val  = q_init;
p.qd_des.val  = zeros(pb.model.NB,1);

p.r_des.val = zeros(3,1);
p.rd_des.val = zeros(3,1);
p.rdd_des.val = zeros(3,1);
%%


% append to pb structure
pb.p = p;

