clear
close all

%% GENERAL USER SET PARAMETERS
path = '../casadi'; %CASADI PATH
verbose = false; %plot for all gear ratios
save_plot = true;
plot_best_gear = false;

%% Dependencies
restoredefaultpath               % "clean slate" for your matlab path
addpath(genpath(path)) % make sure you have added your OS-specific casadi folder to MATLAB-Optimization
import casadi.*

%% Create struct of important params
% motor related
params.motor_base_torque = 2.75; %Nm based on mini cheetah motor, saturation torque
params.motor_base_free_speed = 190; %rads per second mini cheetah motor
params.motor_torque_intercept = 3.677777777777778; %y-intercept of the power line Nm

% mini cheetah parameters
params.body_width = 0.45;
params.body_mass = 3.5; % mass of body in kg

% trajectory related
step_scaling = 3;
params.N  = 50*step_scaling;   % number of control intervals
params.dt = 0.025/step_scaling; % dynamics dt
params.T  = params.N*params.dt; % duration of stance phase
params.g  = -9.81; 

% functions
params.projectile_motion = @(t,y,v,a) y + v*t + 0.5*a*t^2; % anonymous function for projectile motion

%% Loop over gear ratios
% if weight included, gear ratio cannot start at 1 (impossible to solve)
gear_ratios = 5:0.5:40;

max_heights = zeros(length(gear_ratios), 1);
best_gear = "";
best_height = 0;
best_gear_params = params;

for i = 1:length(gear_ratios)
    % set gear ratio in params
    params.gear_ratio = gear_ratios(i);
    
    %actuator torque-speed approx
    params.p1 = [0,params.motor_torque_intercept*params.gear_ratio];
    params.p2 = [params.motor_base_free_speed/params.gear_ratio, 0];
    params.m_ts = (params.p2(2)-params.p1(2)) / (params.p2(1)-params.p1(1));
    params.b_ts = params.motor_torque_intercept*params.gear_ratio;

    % derive dynamics
    [kinematics,dynamics] = derive_leg(params.gear_ratio, params.body_mass/2); 
    
    % formulate optimization via trapezoidal collocation
    [X, U, sol] = optimize_trajectory(kinematics, dynamics, params, path);
     
    % simulate forward the dynamics
    dt = params.dt;
    t = 0:dt:(params.N*dt);
    z = sol.value(X);
    
    %simulate forward projectile motion to find max height
    terminal_COM_sol = kinematics.COM(z(:,end));
    yi = 0;
    vi    = terminal_COM_sol(4);
    t_peak  = -vi / params.g;  
    zf = z(:,end);
    thetaf = zf(2);
   
    max_heights(i) = full(params.projectile_motion(full(t_peak),yi,vi,params.g));
    if max_heights(i) > best_height
        best_gear = params.gear_ratio;
        best_gear_params = params;
        best_gear_params.t = t;
        best_gear_params.t_peak = t_peak;
        best_gear_params.yi = yi;
        best_gear_params.vi = vi;
        best_gear_params.z = z;
        best_gear_params.thetaf = thetaf;
        best_gear_params.kinematics = kinematics;
        best_gear_params.U = U;
        best_gear_params.sol = sol;
    end
    
    % optionally create graphics of leg jumping and torque plot
   if verbose
       params.t = t;
       params.t_peak = t_peak;
       params.yi = yi;
       params.vi = vi;
       params.z = z;
       params.thetaf = thetaf;
       create_graphics(kinematics, params, U, sol, path)
   end
end

% create gear ratio vs max height plot
fig = figure;
hold on
scatter(gear_ratios, max_heights, 'filled');
title("Peak Jump Height vs Gear Ratio for 1D Hopper")
xlabel("Gear Ratio");
ylabel("Peak Height (m)");
grid on
if save_plot
    first = gear_ratios(1)
    step =  gear_ratios(2)-gear_ratios(1)
    last = step*length(gear_ratios)+gear_ratios(1)-1
    filename = "/Graphics/heightVsGearWithWeight"+ first + ":" + step + ":"+ last +".png";
    saveas(fig, pwd + filename);
end

if plot_best_gear
    create_graphics(best_gear_params.kinematics, best_gear_params. best_gear_params.U, best_gear_params.sol, path);
end

    
% Fiona's changes to make:
% fix coloring
% 3D plot, gear_ratio, mass (0:0.1:7), max height
    % actually just mass vs optimal gear ratio
    
    
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

   
