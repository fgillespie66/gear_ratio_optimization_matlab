clear
close all

%% Dependencies
restoredefaultpath               % "clean slate" for your matlab path
addpath(genpath('../casadi')) % make sure you have added your OS-specific casadi folder to MATLAB-Optimization
import casadi.*

%% Create struct of important params
% motor related
params.motor_base_torque = 2.75; %Nm based on mini cheetah motor, saturation torque
params.motor_base_free_speed = 190; %rads per second mini cheetah motor
params.motor_torque_intercept = 3.677777777777778; %y-intercept of the power line Nm

% trajectory related
params.N  = 50;   % number of control intervals
params.dt = 0.025; % dynamics dt
params.T  = params.N*params.dt; % duration of stance phase
params.g  = -9.81; 

% functions
params.projectile_motion = @(t,y,v,a) y + v*t + 0.5*a*t^2; % anonymous function for projectile motion

%% Loop over gear ratios
gear_ratios = [6, 8];


for gear_ratio = gear_ratios
    % set gear ratio in params
    params.gear_ratio = gear_ratio;
    
    % derive dynamics
    [kinematics,dynamics] = derive_leg(gear_ratio); 
    
    % formulate optimization via trapezoidal collocation
     [X, U, sol] = optimize_trajectory(kinematics, dynamics, params, gear_ratio);
     
     % simulate forward the dynamics
    create_graphics(kinematics, params, X, U, sol)
end
    
    
    
    
    
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
