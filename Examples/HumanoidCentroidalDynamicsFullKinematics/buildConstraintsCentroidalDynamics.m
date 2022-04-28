function pb = buildConstraintsCentroidalDynamics(pb)

disp_box('Building constraints');
tic;

%unpack decision variables
var = pb.var;
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
opt_p = pb.p;
m = pb.model.mass;
nDOF = pb.model.NB;
cs = opt_p.cs.sym;
mu = pb.model.mu;
model=pb.model;

opti = pb.opti;
N  = pb.N;
dt = pb.p.dt.sym;

fx_sel = [1 4];
fy_sel = [2 5];
fz_sel = [3 6];
%%
% Loop through timesteps, applying constraints
for k = 1:N
    disp([num2str(k) ' of ' num2str(pb.N)]);
    
    % the 'k' suffix indicates the value of the variable at the current
    % timestep
    qk = q(:,k);
    qdotk = qdot(:,k);
    rpyk = q(4:6,k);
    rk = r(:,k);
    rdk = rdot(:,k);
    rddk = rddot(:,k);
    ftoek = ftoe(:,k);
    ctoek = ctoe(:,k);
    fheelk = fheel(:,k);
    cheelk = cheel(:,k);
    hdk = hdot(:,k);
    hk = h(:,k);
    csk = cs(:,k);
    
    com_k = pb.fun.get_COM(qk);
    
    % Centroid momentum
    Ag_k = pb.fun.get_CMM(qk);
    
    % Dynamics
    opti.subject_to(m*rddk == sum(reshape(ftoek,3,model.NLEGS),2)+...
        sum(reshape(fheelk,3,model.NLEGS),2)+m*model.gravity'); % Eq (7a)
    
    hd_temp = pb.fun.get_hdot([ctoek;cheelk],[ftoek;fheelk],rk);
    
    
    opti.subject_to(hdk == hd_temp); % Eq (7b)
    opti.subject_to(hk == Ag_k(1:3,:)*qdotk); % Eq (7c)
    
    
    % Integrate dynamics
    R_body_to_world = rpyToRotMat(qk(4:6));
    
    
    if (k < N)
        
        % assume qdot piecewise linear. Then q is piecewise quadratic
        acc_fb_lin = (qdot(4:6,k+1) - qdotk(4:6))/dt(k);
        opti.subject_to(q(1:3,k+1) - qk(1:3) == qdotk(4:6) * dt(k) + 0.5*acc_fb_lin*dt(k)^2); % Eq (7d)
        
        acc_fb_ang =  (qdot(1:3,k+1)-qdotk(1:3))/dt(k);
        opti.subject_to(q(4:6,k+1) - qk(4:6) == Binv(rpyk)*(R_body_to_world*qdotk(1:3)*dt(k) +...
            R_body_to_world*acc_fb_ang* dt(k)^2)); % Eq (7d)
        
        acc_joints = (qdot(7:end,k+1)-qdotk(7:end))/dt(k);
        opti.subject_to(q(7:end,k+1) - qk(7:end) == qdotk(7:end)*dt(k) + 0.5*acc_joints*dt(k)^2); % Eq (7d)
        
        
        % assume rdot piecewise linear. Then r is piecewise quadratic and
        % rdd is zero-order hold
        % forward euler for the dynamics
        opti.subject_to(h(:,k+1) - hk == hdk * dt(k)); % Eq (7e)
        opti.subject_to(r(:,k+1) - rk == rdk*dt(k) + 0.5*rddk*dt(k)^2); % Eq (7f)
        opti.subject_to(rdot(:,k+1) - rdk == rddk * dt(k)); % Eq (7g)
    end
    
    % Kinematic Constraints
    opti.subject_to(rk == qk(1:3) + R_body_to_world*(com_k)); % Eq (7h)
    
    
    
    foot_posk = get_gc_position(model, qk,[pb.model.gc_indx.toe pb.model.gc_indx.heel]);
    opti.subject_to(ctoek == foot_posk(1:6)); % Eq (7i)
    opti.subject_to(cheelk == foot_posk(7:12)); % Eq (7i)
    
    
    % Contact Constraint
    addContactConstraint(pb,k)
    
    % Non-negative GRF, Eq (8c)
    opti.subject_to(ftoek(fz_sel) >= zeros(2,1));
    opti.subject_to(fheelk(fz_sel) >= zeros(2,1));
    
    % feet always above the ground
    ground_height_tol = 0.004; % can be up to 2mm below ground to accomodate simulator
    opti.subject_to(ctoek(fz_sel) >=  (-ground_height_tol).*ones(2,1));
    opti.subject_to(cheelk(fz_sel) >= (-ground_height_tol).*ones(2,1));
    
    % State & velocity bounds, Eq (7k)
    % activate/deactivate arm according to parameter flag value
    
    if k==1
        
        opti.subject_to(qk == pb.p.q_init.sym);
        opti.subject_to(qdotk == opt_p.qd_init.sym);
    else
        opti.subject_to(qk <= opt_p.q_max.sym);
        opti.subject_to(qk >= opt_p.q_min.sym);
        opti.subject_to(qdotk <= opt_p.qd_max.sym);
        opti.subject_to(qdotk >= opt_p.qd_min.sym);
    end
    
    % Friction Constraints, Eq (7k)
    opti.subject_to(ftoek(fx_sel) <= 0.71*mu*ftoek(fz_sel));
    opti.subject_to(ftoek(fx_sel) >= -0.71*mu*ftoek(fz_sel));
    opti.subject_to(ftoek(fy_sel) <= 0.71*mu*ftoek(fz_sel));
    opti.subject_to(ftoek(fy_sel) >= -0.71*mu*ftoek(fz_sel));
    
    opti.subject_to(fheelk(fx_sel) <= 0.71*mu*fheelk(fz_sel));
    opti.subject_to(fheelk(fx_sel) >= -0.71*mu*fheelk(fz_sel));
    opti.subject_to(fheelk(fy_sel) <= 0.71*mu*fheelk(fz_sel));
    opti.subject_to(fheelk(fy_sel) >= -0.71*mu*fheelk(fz_sel));
    
    % Constrain flight GRF to zero when not in contact, Eq (8a)
    opti.subject_to(ftoek(fz_sel) <= csk([2 4]).*opt_p.fz_max.sym);
    opti.subject_to(fheelk(fz_sel) <= csk([1 3]).*opt_p.fz_max.sym);
    
end

elapsed_time = toc;
fprintf('built general constraints in %0.2f s\n',elapsed_time)
