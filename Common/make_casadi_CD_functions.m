function [get_COM,get_CMM,get_hdot,get_h,get_r,get_rdd,get_rd,...
    get_h_arms,get_h_body_legs,getCRBI_MX] = make_casadi_CD_functions(model,use_jit)
opts=struct();

if use_jit
    opts = append_jit_options(opts,'-O3');
end

m = model.mass;
q_MX = casadi.MX.sym('q_MX',model.NB);
R_body_to_world_MX = rpyToRotMat(q_MX(4:6));
% COM and CMM

getMassMatrix_MX = make_mass_matrix_MX(model);
getCRBI_MX = make_CRBI_MX(model,getMassMatrix_MX);
get_COM = make_COM_MX(model,getCRBI_MX);
get_CMM = make_cmm_MX(model,getMassMatrix_MX,getCRBI_MX,get_COM,opts);

com_MX  = get_COM(q_MX);
compos_MX = q_MX(1:3) + R_body_to_world_MX*(com_MX); % absolute, world frame

Ag_MX = get_CMM(q_MX);

% Angular Momentum Dynamics
c_MX = casadi.MX.sym('c_MX',3*model.N_GC);
f_MX = casadi.MX.sym('f_MX',3*model.N_GC);
r_MX = casadi.MX.sym('r_MX',3);

hdot_MX = casadi.MX.zeros();
for leg = 1:model.N_GC
    hdot_MX = hdot_MX+...
        cross(c_MX(3*(leg-1)+1:3*(leg-1)+3)-r_MX,f_MX(3*(leg-1)+1:3*(leg-1)+3));
end

qd_MX = casadi.MX.sym('qd_MX',model.NB);
h_MX = Ag_MX(1:3,:)*qd_MX;

h_arms_MX = Ag_MX(1:3,model.indx.arm)*qd_MX(model.indx.arm);
h_body_legs_MX = Ag_MX(1:3,[1:6 model.indx.leg])*qd_MX([1:6 model.indx.leg]);


rd_MX = Ag_MX(4:6,:)*qd_MX/m;

rdd_MX = (sum(reshape(f_MX(1:6),3,model.NLEGS),2)+...
    sum(reshape(f_MX(7:12),3,model.NLEGS),2)+m*model.gravity')/m;



% get_COM = casadi.Function('get_COM_MX',{q_MX},{com_MX}); % com rel
% get_CMM = casadi.Function('get_CMM_MX',{q_MX},{Ag_MX});
get_hdot= casadi.Function('get_hdot_MX',{c_MX,f_MX,r_MX},{hdot_MX},opts);
get_h   = casadi.Function('get_h_MX',{q_MX,qd_MX},{h_MX},opts); % angular momentum
get_h_arms   = casadi.Function('get_h_arms_MX',{q_MX,qd_MX},{h_arms_MX},opts); % arm angular momentum
get_h_body_legs   = casadi.Function('get_h_body_legs_MX',{q_MX,qd_MX},{h_body_legs_MX},opts); % torso+legs angular momentum


get_r   = casadi.Function('get_r_MX',{q_MX},{compos_MX}); % abs com position in world frame
get_rd  = casadi.Function('get_rd_MX',{q_MX,qd_MX},{rd_MX}); % abs com vel in world frame
get_rdd = casadi.Function('get_rdd_MX',{f_MX},{rdd_MX}); % abs com accel in world frame

end

function getMassMatrix_MX = make_mass_matrix_MX(model)
    q_MX = casadi.MX.sym('q_MX',model.NB);
    [H_MX,~] = get_mass_matrix(model, q_MX);
    getMassMatrix_MX = casadi.Function('getH',{q_MX},{H_MX});
end

function getCRBI_MX = make_CRBI_MX(model,getMassMatrix_MX)
    q_MX2 = casadi.MX.sym('q_MX2',model.NB);
    U = [eye(6) zeros(6,model.NB-6)];
    H_MX2 = getMassMatrix_MX(q_MX2);
    Ic = U*H_MX2*U';
    getCRBI_MX = casadi.Function('getIC',{q_MX2},{Ic});
end

function get_COM_MX = make_COM_MX(model,getCRBI_MX)
    m = model.mass;

    q_MX = casadi.MX.sym('q_MX',model.NB);
    Ic_MX = getCRBI_MX(q_MX);
    com_MX = (1/m).*[Ic_MX(3,5);Ic_MX(1,6);Ic_MX(2,4)]; %rel
    get_COM_MX = casadi.Function('get_COM_MX',{q_MX},{com_MX}); % com rel

end

function get_CMM = make_cmm_MX(model,getMassMatrix_MX,getCRBI_MX,get_COM_MX,opts)
    q_MX = casadi.MX.sym('q_MX',model.NB);
    Ic_MX = getCRBI_MX(q_MX);
    com_MX = get_COM_MX(q_MX);
    Ag_MX = get_centroid_momentum_matrix(getMassMatrix_MX(q_MX),Ic_MX,com_MX,q_MX(4:6));
    get_CMM = casadi.Function('get_CMM_MX',{q_MX},{Ag_MX},opts);
    
end


