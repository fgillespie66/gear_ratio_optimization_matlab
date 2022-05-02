function [X, U, sol] = optimize_trajectory(kinematics, dynamics, params, gear_ratio, path)

%% Dependencies
restoredefaultpath               % "clean slate" for your matlab path
addpath(genpath(path)) % make sure you have added your OS-specific casadi folder to MATLAB-Optimization
import casadi.*

%% Formulate Optimization
% via trapezoidal Collocation

opti = casadi.Opti(); % Optimization problem

% getting needed constants
N  = params.N;   % number of control intervals
dt = params.dt; % dynamics dt
T  = params.T; % duration of stance phase
g  = params.g;
projectile_motion = params.projectile_motion

motor_base_torque = params.motor_base_torque
motor_base_free_speed = params.motor_base_free_speed
motor_torque_intercept = params.motor_torque_intercept

p1 = params.p1;
p2 = params.p2;
m_ts = params.m_ts;
b_ts = params.b_ts;

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