clear
close all

%% GENERAL USER SET PARAMETERS
path = '../casadi_mac'; %CASADI PATH
verbose = true; %plot for all gear ratios

%% Dependencies
restoredefaultpath               % "clean slate" for your matlab path
addpath(genpath(path)) % make sure you have added your OS-specific casadi folder to MATLAB-Optimization
import casadi.*

%% Create struct of important params
% motor related
params.motor_base_torque = 2.75; %Nm based on mini cheetah motor, saturation torque
params.motor_base_free_speed = 190; %rads per second mini cheetah motor
params.motor_torque_intercept = 3.677777777777778; %y-intercept of the power line Nm

% trajectory related
step_scaling = 3;
params.N  = 50*step_scaling;   % number of control intervals
params.dt = 0.025/step_scaling; % dynamics dt
params.T  = params.N*params.dt; % duration of stance phase
params.g  = -9.81; 

% functions
params.projectile_motion = @(t,y,v,a) y + v*t + 0.5*a*t^2; % anonymous function for projectile motion

%% Loop over gear ratios
gear_ratios = [6, 8];


for gear_ratio = gear_ratios
    % set gear ratio in params
    params.gear_ratio = gear_ratio;
    
    %actuator torque-speed approx
    params.p1 = [0,params.motor_torque_intercept*gear_ratio];
    params.p2 = [params.motor_base_free_speed/gear_ratio, 0];
    params.m_ts = (params.p2(2)-params.p1(2)) / (params.p2(1)-params.p1(1));
    params.b_ts = params.motor_torque_intercept*gear_ratio;

    % derive dynamics
    [kinematics,dynamics] = derive_leg(gear_ratio); 
    
    % formulate optimization via trapezoidal collocation
     [X, U, sol] = optimize_trajectory(kinematics, dynamics, params, gear_ratio, path);
     
     % simulate forward the dynamics
     if verbose
        create_graphics(kinematics, params, X, U, sol, path)
     end
end
    
    
    
 %% TESTING CONSTRAINTS GRAPHICALLY
    %figure
    %hold on
    %plot([0,motor_base_free_speed/gear_ratio],[motor_torque_intercept*gear_ratio,0]); %draw the power limit line
    %fplot(@(x) m_ts * x + b_ts);
    %legend('one', 'two');
%CORRECTION FACTOR TELLS US IT'S JUST THE y-INTERCEPT THAT'S WRONG











%% NOTES FOR OTHER THINGS


%PLAN
%first expand this to simulate the whole leg forward 
%second expand this to a more complex leg model
%third expand this to a mini cheetah 
%CITE the interias paper, yanran's paper, RUS's work, MATT's work, anything
%else we referenced including the Matt Kelly Optimization tutorial, matt
%and charles' PDF tutorial

%remember, feedback control is REQUIRED just to take out errors in your
%modeling!!! 

%can we use a MATLAB function as a constraint? if we write three
%constraints 



%SOME OTHER NOTES
%add the line as a COST? 
%we learned that the constraint doesn't actually work in this case 
% might need a quper quadratic cost ?? 


%print the constraint violations at each time step 

%TRY DECREASING THE TIMESTEP??? 
%problem gets worse @ higher gear ratios, 
%decreasing timestep significantly helped 

%% OK HERE's HOW WE FIXED IT 

%look @ our integration scheme 
% V_(k+1) = Vk + Uk*dt

%so the arrays look like this !!! 
% [v1, v2, v2, v4 ... vN]
% [u1, u2, u3, ... uN-1] 

%the torques go with the FIRST time step because they get APPLIED to
%calculate the next time step 

%which also asks the question should we constrain Uk to Vk or Vk+1 ??? 

   
