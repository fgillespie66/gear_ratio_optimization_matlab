function [] = create_graphics(kinematics, params, X, U, sol)
%% Dependencies
restoredefaultpath               % "clean slate" for your matlab path
addpath(genpath('../casadi')) % make sure you have added your OS-specific casadi folder to MATLAB-Optimization
import casadi.*

%% Get params
N  = params.N;   % number of control intervals
dt = params.dt; % dynamics dt
T  = params.T; % duration of stance phase
g  = params.g;
projectile_motion = params.projectile_motion

motor_base_torque = params.motor_base_torque
motor_base_free_speed = params.motor_base_free_speed
motor_torque_intercept = params.motor_torque_intercept
gear_ratio = params.gear_ratio


%% Simulate Forward the Dynamics

% ---- post-processing        ------
t = 0:dt:(N*dt);
z = sol.value(X);

%---- simulating forward the dyanmics ----
terminal_COM_sol = kinematics.COM(z(:,end));
yi = 0;
vi    = terminal_COM_sol(4);
t_peak   = -vi / g;  
zf = z(:,end);
thetaf = zf(2);

%simulate forward projectile motion 
t2 = dt:dt:full(t_peak);
ys = [];
ts = [];
gs = [];
zeros = [];
for time = t2
    yc = projectile_motion(time,yi,vi,g); %use this with yi = 0 but vi is the same! 
    ys = [ys, full(yc)];
    ts = [ts, thetaf];
    gs = [gs, g];
    zeros = [zeros, 0.0];
end

%% Animate the Solution

%prep the arrays for plotting
t2 = t2+N*dt;
zs = [ys; ts; gs; zeros];
t_fs = [t, t2];
z_fs = [z, zs];

%animate the solution
figure;
animate_simple(t_fs,z_fs,kinematics,1);

%% Plot Actuation Efforts + True Torque Speed Curve

%populate the torque curve and the speed curve 
%motor_torques = abs(sol.value(U)); %take the absolute value of torques
motor_torques = sol.value(U);
motor_velocities = [];
for x = z
    %velocity = abs(x(end)); %take the absolute value of the actuation commands
    velocity = x(end);
    motor_velocities = [motor_velocities, velocity];
end

%remove the initial configuration of the leg from state vector
motor_velocities = motor_velocities(2:end);

%formatting the figure
sz = 25;
c = linspace(1,10,length(motor_velocities));
%plot
figure;
hold on
scatter(motor_velocities, motor_torques, sz, c, 'filled');
plot([0,motor_base_free_speed/gear_ratio],[motor_base_torque*gear_ratio,motor_base_torque*gear_ratio]); %draw the torque saturation limit line
plot([motor_base_free_speed/gear_ratio,motor_base_free_speed/gear_ratio],[0,motor_base_torque*gear_ratio]); %draw the free speed
plot([0,motor_base_free_speed/gear_ratio],[motor_torque_intercept*gear_ratio,0]); %draw the power limit line
xlabel("Motor Velocity (rads/s)");
ylabel("Motor Torques (Nm)");
title("Actuation Command Overlayed on TS Curve during Stance");
legend('Trajectory Start', 'Saturation Torque', 'Free Speed', 'Physical Limit');
grid on
