clear
close all

%% Dependencies
restoredefaultpath               % "clean slate" for your matlab path
addpath(genpath('../casadi_mac')) % make sure you have added your OS-specific casadi folder to MATLAB-Optimization
import casadi.*

%% Derive dynamics

%actuator base parameters
gear_ratio = 12;
motor_base_torque = 2.75; %Nm based on mini cheetah motor, saturation torque
motor_base_free_speed = 190; %rads per second mini cheetah motor
motor_torque_intercept = 3.677777777777778; %y-intercept of the power line Nm

%mini cheetah parameters
body_width = 0.45;
body_mass = 3.5; %mass of body in kg 

%actuator torque-speed approx
p1 = [0,motor_torque_intercept*gear_ratio];
p2 = [motor_base_free_speed/gear_ratio, 0];
m_ts = (p2(2)-p1(2)) / (p2(1)-p1(1));
b_ts = motor_torque_intercept*gear_ratio;

%derive dynamics
[kinematics,dynamics] = derive_leg(gear_ratio, body_mass/2); 

%% Formulate Optimization
% via trapezoidal Collocation

opti = casadi.Opti(); % Optimization problem

step_scaling = 5;
N  = 50*step_scaling;   % number of control intervals
dt = 0.025/step_scaling; % dynamics dt
T  = N*dt; % duration of stance phase

% ---- decision variables ---------
X = opti.variable(4,N+1); % [y, theta, dy, dtheta] NOTE: this includes initial conditions ! 
U = opti.variable(1,N);   % hip torque (could alternatively parameterize U by a spline
F = opti.variable(1,N);   % vertical reaction force

% ---- objective          ---------
g = -9.81;
terminal_COM          = kinematics.COM(X(:,end)); 
terminal_com_y_height = terminal_COM(2);
terminal_com_y_vel    = terminal_COM(4);

time_to_peak   = -terminal_com_y_vel / g;           % find time when
projectile_motion = @(t,y,v,a) y + v*t + 0.5*a*t^2; % anonymous function for projectile motion
max_com_height = projectile_motion(time_to_peak,terminal_com_y_height,terminal_com_y_vel,g);

opti.minimize(-max_com_height); % maximize peak com height
%opti.minimize(-max_com_height-0.1*terminal_com_y_vel); % maximize peak com height

% ---- dynamic constraints --------
for k=1:N % loop over control intervals
    Xk  = X(:,k); 
    Xk1 = X(:,k+1);
    Uk  = U(:,k);
    Fk  = F(:,k);
    Ak1 = dynamics.A(Xk1);
    bk1 = dynamics.b(Xk1, Uk, Fk);
   opti.subject_to( Xk1(1:2) - Xk(1:2) == dt*Xk1(3:4) ) % Euler integration - position
   opti.subject_to( Ak1*(Xk1(3:4)-Xk(3:4))  == dt*bk1 ) % Euler integration - velocity
end

% ---- path constraints -----------
for k=1:N % loop over control intervals
    Xk  = X(:,k); 
    Xk1 = X(:,k+1);
    Uk  = U(:,k);
    Fk  = F(:,k);
    Vm = Xk(end);
    opti.subject_to( Xk1(1) == 0 )      % Foot in place
    opti.subject_to( Fk >= 0 )          % Unilateral force constraint (no pulling the ground)
    opti.subject_to( -motor_base_torque*gear_ratio <= Uk <= motor_base_torque*gear_ratio)   % Control limits INCLUDE GEAR RATIO
    opti.subject_to( -motor_base_free_speed/gear_ratio <= Vm <= motor_base_free_speed/gear_ratio)

    %torque speed curve line torque <= m * velocity + b
    opti.subject_to( (Uk - m_ts * Vm) <= (b_ts) )

    % FROM MATT: TRY TO FORMAT ALL OF THESE CONSTRAINTS AS A SINGLE Ax <= b
    % CONSTRAINT!!! THAT MIGHT MAKE THE OPTIMIZER HAPPIER 

    %add in a limit for the motor velocity (i.e. we cap theta_dot)! 
    opti.subject_to( 0 <= Xk1(2) <= pi/2.2 )% Knee above ground, avoid singularity
end

% ---- boundary conditions --------
opti.subject_to( X(:,1) == [0;deg2rad(10);0;0] );

% ---- initial values for solver ---
%opti.set_initial(X, );
%opti.set_initial(U, );
%opti.set_initial(F, );

% ---- solve NLP              ------
p_opts = struct('expand',true); % expand to casadi variables to SX (10x speedup)
opti.solver('ipopt',p_opts);    % set numerical backend
sol = opti.solve();             % actual solve

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
MAX_HEIGHT = projectile_motion(full(t_peak),yi,vi,g)
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
figure;
animate_simple(t_fs,z_fs,kinematics,1, gear_ratio, body_width);

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

%formatting the figure
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
title("Actuation Command Overlayed on TS Curve during Stance");
legend('Trajectory Start', 'Saturation Torque', 'Free Speed', 'Physical Limit');
grid on










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

