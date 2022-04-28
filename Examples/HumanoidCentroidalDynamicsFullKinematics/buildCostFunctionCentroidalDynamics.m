function pb  = buildCostFunctionCentroidalDynamics(pb)

%unpack decision variables
var = pb.var;
X=var.X;
q = var.q;
qdot = var.qdot;
r = var.r;
rdot = var.rdot;
rddot = var.rddot;
ctoe = var.ctoe;
cheel  = var.cheel;
ftoe = var.ftoe;
fheel=var.fheel;
h = var.h;
hdot = var.hdot;

% unpack other stuff
N = pb.N;
opt_p = pb.p;
dt = opt_p.dt.sym;
stance_time = pb.stance_time;

q_des = pb.p.q_des.sym;
qd_des = pb.p.qd_des.sym;
r_des = pb.p.r_des.sym;
rd_des = pb.p.r_des.sym;
rdd_des = pb.p.r_des.sym;

g = pb.model.g;


cost = casadi.MX.zeros(1,1);
disp_box('building general cost')
% Running cost
for k = 1:N
    
    qk = q(:,k);
    qdotk = qdot(:,k);
    rk = r(:,k);
    rdk = rdot(:,k);
    rddk = rddot(:,k);
    ftoek = ftoe(:,k);
    fheelk = fheel(:,k);
    ctoek = ctoe(:,k);
    cheelk = cheel(:,k);
    hk = h(:,k);
    hdotk = hdot(:,k);
    p_gaitk = pb.p.cs.sym(:,k);
    
    % setup desired command in cost function
    R_yaw_frame_to_world = rpyToRotMat([0;0;qk(6)]);%
    desired_xy_vel_world_frame = R_yaw_frame_to_world(1:2,1:2)*opt_p.desired_state_cmd.sym(2:3); % rotate desired xy velocity in
       
    qk_des = [q_des(1:2);opt_p.desired_state_cmd.sym(4);q_des(4:24)];
    
    qdk_des = [qd_des(1:2);opt_p.desired_state_cmd.sym(1);desired_xy_vel_world_frame;qd_des(6:24)];
    
   
    cref = repmat(qk(1:3),pb.model.N_GC,1) + opt_p.pf_ref.sym;

    href  = [0;0;0];
    hdref = [0;0;0];

    

    fref = [R_yaw_frame_to_world*[0;0;pb.model.mass*g*pb.Tp/pb.stance_time/4].*p_gaitk(2);...
        R_yaw_frame_to_world*[0;0;pb.model.mass*g*pb.Tp/pb.stance_time/4].*p_gaitk(4);...
        R_yaw_frame_to_world*[0;0;pb.model.mass*g*pb.Tp/pb.stance_time/4].*p_gaitk(1);...
        R_yaw_frame_to_world*[0;0;pb.model.mass*g*pb.Tp/pb.stance_time/4].*p_gaitk(3)];

    
    Xref = [qk_des;qdk_des;r_des;rd_des;rdd_des;cref;fref;href;hdref];
    QX = [opt_p.Qq.sym;opt_p.Qv.sym;opt_p.Qr.sym;opt_p.Qc.sym;opt_p.Qf.sym;opt_p.Qh.sym];

    Xerr = X(:,k) - Xref;
    
    if k < N
        cost_new = (Xerr'*diag(QX)*Xerr)*dt(k);
    else
        cost_new = (Xerr'*diag(QX)*Xerr)*dt(end); % patch for terminal cost for now
    end
    
    
    cost = cost + cost_new;
end


% raibert cost
% disp_box('building raibert cost')
% raibert_cost = getRaibertCost(pb,phase);
% cost = cost + raibert_cost;

%%
% attach final cost to pb structure for output
pb.cost = cost;
pb.opti.minimize(cost);