function [kinematics,dynamics] = derive_leg_2DOF(gear_ratio, payload_mass) 
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
y    = casadi.SX.sym('y',1,1);    % vertical position of CoM
th1   = casadi.SX.sym('th',1,1);   % hip angle (from vertical)
th2   = casadi.SX.sym('th',1,1);   % knee (from hip)
dy   = casadi.SX.sym('dy',1,1);   % vertical velocity
dth1  = casadi.SX.sym('dth',1,1);  
dth2  = casadi.SX.sym('dth',1,1); 
ddy  = casadi.SX.sym('ddy',1,1);  % vertical acceleration
ddth1 = casadi.SX.sym('ddth',1,1); % shin angular acceleration
ddth2 = casadi.SX.sym('ddth',1,1); % shin angular acceleration
Fy   = casadi.SX.sym('Fy',1,1);   % constraint force
tau1  = casadi.SX.sym('tau',1,1);  % hip torque
tau2  = casadi.SX.sym('tau',1,1);  % hip torque

%MINI CHEETAH PARAMETERS 
l = 0.209;
c1 = 0.098;
c2 = 0.209/2;
m1 = 0.092;   % link 1 mass
m2 = 0.06;   % link 2 mass
mh = 0.440+payload_mass/2;   % hip mass
I1 = 0.05;  % link 1 inertia
I2 = 0.05;  % link 2 inertia
Ir = 33*10^(-6);  % rotor inertia kgm^2
N  = gear_ratio;     % gear ratio
g = 9.81;

% Group terms for later use
q   = [y; th1; th2];      % generalized coordinates
dq  = [dy; dth1; dth2];    % first time derivatives
ddq = [ddy; ddth1; ddth2];  % second time derivatives
u   = [tau1, tau2];          % control forces and moments
Fc   = Fy;          % constraint forces and moments

%%% Calculate important vectors and their time derivatives.

% Define fundamental unit vectors.  The first element should be the
% horizontal (+x cartesian) component, the second should be the vertical (+y
% cartesian) component, and the third should be right-handed orthogonal.
ihat = [1; 0; 0];
jhat = [0; 1; 0];
khat = cross(ihat,jhat);

% Define other unit vectors for use in defining other vectors.
er1hat =  cos(th1)*ihat - sin(th1) * jhat;
er2hat =  cos(th2)*ihat - sin(th2) * jhat;

% A handy anonymous function for taking first and second time derivatives
% of vectors using the chain rule.  See Lecture 6 for more information. 
ddt = @(r) jacobian(r,[q;dq])*[dq;ddq]; 

% Define vectors to key points.
d = l*sin(th1)+l*sin(th2); %height of leg at any given time
rh   = (d+y)*jhat;         % hip
rcm1 = rh+c1*er1hat;   % COM of link 1
rk   = rh+l*er1hat;    % knee
rcm2 = rk + c2*er2hat; % COM of link 2
rf   = y*jhat+(l*cos(th1)+l*cos(th2))*ihat;  % foot
%rf   = y*jhat;
keypoints = [rf rk rh];

% Take time derivatives of vectors as required for kinetic energy terms.
drcm1 = ddt(rcm1);
drcm2 = ddt(rcm2);
drf = ddt(rf);
drh  = ddt(rh);

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
T1  = (1/2)*m1*dot(drcm1, drcm1) + (1/2)* I1 * dth1^2;   % link 1 KE
T2  = (1/2)*m2*dot(drcm2, drcm2) + (1/2)* I2 * dth2^2;   % link 2 KE
T2r = (1/2)*Ir*(dth2 + N*dth2)^2;                          % link 2 rotor inertia
T1r = (1/2)*Ir*(dth1 + N*dth1 + dth2)^2;                   % link 1 rotor inertia include rotor of link 2
Tcom  = (1/2)*mh*dot(drh, drh);                            % hip KE

% Define potential energies. See Lecture 6 formulas for gravitational 
% potential energy of rigid bodies and elastic potential energies of
% energy storage elements.
V1 = m1*g*dot(rcm1, jhat);
V2 = m2*g*dot(rcm2, jhat);
Vh = mh*g*dot(rh, jhat);

% Define contributions to generalized forces.  See Lecture 6 formulas for
% contributions to generalized forces.

% THIS WE NEED TO CHANGE~~~

%QF   = F2Q(Fy*jhat,rf);
%Qtau = M2Q(-tau*khat, -dth*khat);

% Sum kinetic energy terms, potential energy terms, and generalized force
% contributions.
T = T1 + T1r + T2 + T2r + Tcom;
V = V1 + V2 + Vh;
Q = [Fy; tau1; tau2];
%Q = QF + Qtau;

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
