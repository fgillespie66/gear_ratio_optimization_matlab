function addContactConstraint(pb,k)

%unpack decision variables
var = pb.var;
ctoe = var.ctoe;
cheel  = var.cheel;

cs = pb.p.cs.sym;
N  = pb.N;
ctoek = ctoe(:,k);
cheelk = cheel(:,k);
csk = cs(:,k);
opti = pb.opti;

% Contact Constraint
for leg = 1:pb.model.NLEGS
    % if a heel is in contact
    
    opti.subject_to(cs(2*leg-1,k).*cheelk(3*(leg-1)+3) == 0);
    
    if k < N
          stay_on_ground = repmat(cs(2*leg-1,k),2,1);%if next time step also a contact
        opti.subject_to(stay_on_ground.*(cheel(3*(leg-1)+1:3*(leg-1)+2,k+1)-cheelk(3*(leg-1)+1:3*(leg-1)+2)) == 0); % see Gerardo's thesis p.147
    end
    
    % if a toe is in contact
    opti.subject_to(cs(2*leg,k).*ctoek(3*(leg-1)+3) == 0);
    
    if k < N
        stay_on_ground = repmat(cs(2*leg,k),2,1);
        opti.subject_to(stay_on_ground.*(ctoe(3*(leg-1)+1:3*(leg-1)+2,k+1)-ctoek(3*(leg-1)+1:3*(leg-1)+2)) == 0); % see Gerardo's thesis p.147
    end
end