function [] = create_graphics(kinematics, params, U, sol, path)

%% Dependencies
restoredefaultpath               % "clean slate" for your matlab path
addpath(genpath(path)) % make sure you have added your OS-specific casadi folder to MATLAB-Optimization
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

p1 = params.p1;
p2 = params.p2;
m_ts = params.m_ts;
b_ts = params.b_ts;

t = params.t;
z = params.z;
t_peak = params.t_peak;
yi = params.yi;
vi = params.vi;
thetaf = params.thetaf;

%% Simulate Forward the Dynamics
t2 = dt:dt:2*full(t_peak); %simulate to right before impact
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
fig = figure;
animate_simple(t_fs,z_fs,kinematics,1, gear_ratio);


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
motor_velocities = motor_velocities(1:end-1); %DONT REMOVE THE FIRST ONE

%formatting the figure;
sz = 25;
c = linspace(1,10,length(motor_velocities));
%plot
figure;
hold on
scatter(motor_velocities, motor_torques, sz, c, 'filled');
%plot(motor_velocities, motor_torques); MAYBE ADD A LINE HERE FOR
%orientation" also change the thicknesses of everything
plot([0,motor_base_free_speed/gear_ratio],[motor_base_torque*gear_ratio,motor_base_torque*gear_ratio]); %draw the torque saturation limit line
plot([motor_base_free_speed/gear_ratio,motor_base_free_speed/gear_ratio],[0,motor_base_torque*gear_ratio]); %draw the free speed
%plot([0,motor_base_free_speed/gear_ratio],[motor_torque_intercept*gear_ratio,0]); %draw the power limit line
fplot(@(x) m_ts * x + b_ts); %draw the power limit line
xlabel("Motor Velocity (rads/s)");
ylabel("Motor Torques (Nm)");
title("Gear Ratio: " + gear_ratio + " - Stance Actuation Overlayed on TS Curve");
legend('Trajectory Start', 'Saturation Torque', 'Free Speed', 'Physical Limit');
grid on