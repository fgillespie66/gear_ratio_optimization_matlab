function pb = make_parameters(pb)

% unpack structure elements
nDOF = pb.model.NB;
opti = pb.opti;
N = pb.N;

if isfield(pb,'p') % if p field already exists
    p = pb.p; % keep already existing structure parameters
end

% other parameters
% p.initial_state.sym = opti.parameter(2*nDOF,1);
p.desired_state_cmd.sym = opti.parameter(4,1);
p.raibert.sym = opti.parameter(2,1);

% cost function weights (Phase 1)
% store the diagonal of cost matrices here
p.Qq.sym  = opti.parameter(nDOF,1);
p.Qv.sym  = opti.parameter(nDOF,1);
p.Qr.sym  = opti.parameter(9,1);
p.Qh.sym  = opti.parameter(6,1);
p.Qc.sym = opti.parameter(12,1);
p.Qf.sym  = opti.parameter(12,1);

p.Q_step_reg.sym = opti.parameter(4,1);

% timing parameters
p.dt.sym = opti.parameter(1,pb.N-1); % leave as constant scalar for now
p.cs.sym = opti.parameter(pb.model.N_GC,pb.N);

% state bounds
p.q_min.sym  = opti.parameter(nDOF,1);
p.q_max.sym  = opti.parameter(nDOF,1);
p.qd_min.sym = opti.parameter(nDOF,1);
p.qd_max.sym = opti.parameter(nDOF,1);

p.q_init.sym  = opti.parameter(nDOF,1); % initial state not used yet TODO Use this instead
% p.q_init_max.sym  = opti.parameter(nDOF,1);
p.qd_init.sym = opti.parameter(nDOF,1);
% p.qd_init_max.sym = opti.parameter(nDOF,1);

p.q_final_min.sym  = opti.parameter(nDOF,1);
p.q_final_max.sym   = opti.parameter(nDOF,1);
p.qd_final_min.sym = opti.parameter(nDOF,1);
p.qd_final_max.sym = opti.parameter(nDOF,1);
p.pf_ref.sym = opti.parameter(12,1);


% force bounds
p.fz_max.sym = opti.parameter(1,1);
p.fz_min.sym = opti.parameter(1,1);

%reference

p.q_des.sym = opti.parameter(pb.model.NB,1);
p.qd_des.sym = opti.parameter(pb.model.NB,1);

p.r_des.sym = opti.parameter(3,1);
p.rd_des.sym = opti.parameter(3,1);
p.rdd_des.sym = opti.parameter(3,1);
%%
pb.p = p;


end