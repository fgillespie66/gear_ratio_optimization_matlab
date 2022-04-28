function pb = make_decision_variables(pb)
% centroidal dynamics and full kinematics
% Optimization variables
% q       (nDOF) [x,y,z,roll,pitch,yaw,joint pos]
% qdot    (nDOF) [omega (body frame), vBody (body frame), joint vel]
% r       (3)
% rdot    (3)
% rddot   (3)
% ctoe    (NLEGS*3)
% cheel   (NLEGS*3)
% ftoe    (NLEGS*3)
% fheel   (NLEGS*3)
% h       (3)
% hdot    (3)

nDOF = pb.model.NB;

X = pb.opti.variable(2*nDOF+39, pb.N);
var.X=X;
%unwrapped decision variables
var.q       = X(1:nDOF,:);              % fb state + joint state
var.qdot    = X(nDOF+1:2*nDOF,:);       % fb velocity + joint velocity
var.r       = X(2*nDOF+1:2*nDOF+3,:);   % com position
var.rdot    = X(2*nDOF+4:2*nDOF+6,:);   % com velocity
var.rddot   = X(2*nDOF+7:2*nDOF+9,:);   % com acceleration
var.ctoe    = X(2*nDOF+10:2*nDOF+15,:); % toe positions
var.cheel   = X(2*nDOF+16:2*nDOF+21,:); % heel positions
var.ftoe    = X(2*nDOF+22:2*nDOF+27,:); % ground reaction forces
var.fheel   = X(2*nDOF+28:2*nDOF+33,:); % ground reaction forces
var.h       = X(2*nDOF+34:2*nDOF+36,:); % centroid angular momentum
var.hdot    = X(2*nDOF+37:2*nDOF+39,:); % deriv of CAM

pb.var = var;
end