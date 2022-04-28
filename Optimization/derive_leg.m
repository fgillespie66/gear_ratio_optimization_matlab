function [kinematics,dynamics] = derive_leg(gear_ratio) 
% Two link hopper (need to draw a diagram of the leg to reference)
%
% Foot and hip are constrained to a vertical rail
%
% Only hip is actuated
% 
% Returns 2 structures
% @kinematics: contains casadi functions for kinematics
% @dynamics: contains casadi functions for dynamics

% Define variables for generalized coordinates + derivatives and controls
y    = casadi.SX.sym('y',1,1);    % vertical position
th   = casadi.SX.sym('th',1,1);   % shin angle
dy   = casadi.SX.sym('dy',1,1);   % vertical velocity
dth  = casadi.SX.sym('dth',1,1);  % shin angular velocity
ddy  = casadi.SX.sym('ddy',1,1);  % vertical acceleration
ddth = casadi.SX.sym('ddth',1,1); % shin angular acceleration
tau  = casadi.SX.sym('tau',1,1);  % hip torque (knee is not actuated)
Fy   = casadi.SX.sym('Fy',1,1);   % constraint force

% Parameters (original)
% l = 0.5;
% c1 = 0.25;
% c2 = 0.25;
% m1 = 0.1;   % link 1 mass
% m2 = 0.1;   % link 2 mass
% mh = 0.3;   % hip mass
% I1 = 0.05;  % link 1 inertia
% I2 = 0.05;  % link 2 inertia
% Ir = 0.02;  % rotor inertia kgm^2
% N  = gear_ratio;     % gear ratio
% g = 9.81;

%MINI CHEETAH PARAMETERS 
l = 0.209;
c1 = 0.098;
c2 = 0.209/2;
m1 = 0.092;   % link 1 mass
m2 = 0.06;   % link 2 mass
mh = 0.440;   % hip mass
I1 = 0.05;  % link 1 inertia
I2 = 0.05;  % link 2 inertia
Ir = 33*10^(-6);  % rotor inertia kgm^2
N  = gear_ratio;     % gear ratio
g = 9.81;

% Group terms for later use
q   = [y; th];      % generalized coordinates
dq  = [dy; dth];    % first time derivatives
ddq = [ddy; ddth];  % second time derivatives
u   = tau;          % control forces and moments
Fc   = Fy;          % constraint forces and moments
%p   = [l; c1; c2; m1; m2; mh; I1; I2; g];  % parameters

%%% Calculate important vectors and their time derivatives.

% Define fundamental unit vectors.  The first element should be the
% horizontal (+x cartesian) component, the second should be the vertical (+y
% cartesian) component, and the third should be right-handed orthogonal.
ihat = [1; 0; 0];
jhat = [0; 1; 0];
khat = cross(ihat,jhat);

% Define other unit vectors for use in defining other vectors.
er1hat =  cos(th)*ihat + sin(th) * jhat;
er2hat = -cos(th)*ihat + sin(th) * jhat;

% A handy anonymous function for taking first and second time derivatives
% of vectors using the chain rule.  See Lecture 6 for more information. 
ddt = @(r) jacobian(r,[q;dq])*[dq;ddq]; 

% Define vectors to key points.
rf   = y*jhat;         % foot
rcm1 = rf+c1*er1hat;   % COM of link 1
rk   = rf+l*er1hat;    % knee
rcm2 = rk + c2*er2hat; % COM of link 2
rh   = rk + l*er2hat;  % hip
keypoints = [rh rk rf];

% Take time derivatives of vectors as required for kinetic energy terms.
drcm1 = ddt(rcm1);
drcm2 = ddt(rcm2);
drh   = ddt(rh);

%%% Calculate Kinetic Energy, Potential Energy, and Generalized Forces

% F2Q calculates the contribution of a force to all generalized forces
% for forces, F is the force vector and r is the position vector of the 
% point of force application
F2Q = @(F,r) simplify(jacobian(r,q)'*(F)); 

% M2Q calculates the contribution of a moment to all generalized forces
% M is the moment vector and w is the angular velocity vector of the
% body on which the moment acts
M2Q = @(M,w) simplify(jacobian(w,dq)'*(M)); 

% Define kinetic energies. See Lecture 6 formula for kinetic energy
% of a rigid body.
T1  = (1/2)*m1*dot(drcm1, drcm1) + (1/2)* I1 * dth^2;    % link 1 KE
T2  = (1/2)*m2*dot(drcm2, drcm2) + (1/2)* I2 * (-dth)^2; % link 2 KE
T2r = (1/2)*Ir*(dth + N*dth)^2;                          % link 2 rotor inertia
Th  = (1/2)*mh*dot(drh, drh);                            % hip KE

% Define potential energies. See Lecture 6 formulas for gravitational 
% potential energy of rigid bodies and elastic potential energies of
% energy storage elements.
V1 = m1*g*dot(rcm1, jhat);
V2 = m2*g*dot(rcm2, jhat);
Vh = mh*g*dot(rh, jhat);

% Define contributions to generalized forces.  See Lecture 6 formulas for
% contributions to generalized forces.
QF   = F2Q(Fy*jhat,rf);
Qtau = M2Q(-tau*khat, -dth*khat);

% Sum kinetic energy terms, potential energy terms, and generalized force
% contributions.
T = T1 + T2 + T2r + Th;
V = V1 + V2 + Vh;
Q = QF + Qtau;

% Calculate rcm, the location of the center of mass
rcm = (m1*rcm1 + m2*rcm2 + mh*rh)/(m1+m2+mh);

% Assemble C, the set of constraints
C = y;  % When y = 0, the constraint is satisfied because foot is on the ground
dC= ddt(C);

%% All the work is done!  Just turn the crank...
%%% Derive Energy Function and Equations of Motion
E = T+V;                                         % total system energy
L = T-V;                                         % the Lagrangian
eom = ddt(jacobian(L,dq)') - jacobian(L,q)' - Q;  % form the dynamics equations

%size(eom)

%%% Rearrange Equations of Motion. 
A = jacobian(eom,ddq);
b = A*ddq - eom;

%%% Write functions to evaluate dynamics, etc...
z = [q;dq];
dynamics.energy      = casadi.Function('energy',{z},{E});
dynamics.A           = casadi.Function('A',{z},{A});
dynamics.b           = casadi.Function('b',{z u Fc},{b});
kinematics.keypoints = casadi.Function('keypoints',{z},{keypoints});
kinematics.C         = casadi.Function('constraints',{z u},{C});
kinematics.dC        = casadi.Function('dconstraints',{z u},{dC});

% Write a function to evaluate the X and Y coordinates and speeds of the center of mass given the current state and parameters
drcm = ddt(rcm);             % Calculate center of mass velocity vector
COM = [rcm(1:2); drcm(1:2)]; % Concatenate x and y coordinates and speeds of center of mass in array
kinematics.COM = casadi.Function('COM',{z},{COM});
