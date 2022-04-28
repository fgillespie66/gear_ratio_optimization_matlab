clear;close all;clc;
%% Description
% This script implement trajectory optimization based on centroidal dynamics and full kinematics
% Dai, Hongkai, Andres Valenzuela, and Russ Tedrake. “Whole-Body
% Motion Planning with Centroidal Dynamics and Full Kinematics.”
% 2014 IEEE-RAS International Conference on Humanoid Robots
% (November 2014).
% see https://dspace.mit.edu/handle/1721.1/101079
%% Add relevant libraries
% Restore default path
restoredefaultpath
% Libraries for all OS
addpath(genpath('.'));% current dir
addpath(genpath('../../spatial_v2_extended'));
addpath(genpath('../../casadi'));
addpath(genpath('../../FloatingBaseDynamics3D'));
addpath(genpath('../../Robots'));
addpath(genpath('../../Common'));
addpath(genpath('../../Animation'));

% OS Specific Libraries
import casadi.*
%% Build model
disp_box('Building Robot Model');
pb.model = getHumanoidModel('eul');
pb.model.mu = 0.7; % friction coefficient
pb.model.g = abs(pb.model.gravity(3));

[pb.fun.get_COM,pb.fun.get_CMM,pb.fun.get_hdot,pb.fun.get_h,pb.fun.get_r,pb.fun.get_rdd,pb.fun.get_rd,...
    ~,~,~] = make_casadi_CD_functions(pb.model,false);

%% Timing
T = 0.8; % time horizon;
dt = 0.1; % discretization dt
N = floor(T/dt) + 1;

pb.N = N;
pb.stance_time = T/2;
pb.Tp = T; % gait period time
pb.dt = dt;

%% Create decision variables and optimization parameters
pb.opti = casadi.Opti();
pb = make_decision_variables(pb);
pb = make_parameters(pb);
pb = buildCostFunctionCentroidalDynamics(pb);
pb = buildConstraintsCentroidalDynamics(pb);


%% get parameters values and set them
pb = getParameterValuesCD(pb);
pb = setParameterValues(pb);

%% Set initial guess
% TODO should set initial guess for all decision variables. if no guess is
% set, default is 0;

q_x_guess = linspace(pb.p.q_init.val(1),pb.p.q_init.val(1)+T*pb.p.desired_state_cmd.val(2),pb.N);
pb.opti.set_initial(pb.var.q(1,:), q_x_guess)
%% Solve !
p_opts = struct('expand',true); % this speeds up ~x10
s_opts = get_solver_options('ipopt');


pb.opti.solver('ipopt',p_opts,s_opts);

disp_box(['Solving']);
try
    sol = pb.opti.solve_limited();
catch ME
    try
        pb.opti.debug.show_infeasibilities(1e-6)
    catch
        rethrow(ME)
    end
end

%% plot solution
q_star = pb.opti.value(pb.var.q);
t_star = 0:dt:(pb.N-1)*dt;

model = buildShowMotionModel(pb.model, 0);

% Animation
showmotion(model,t_star,q_star)