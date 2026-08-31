function par = transformParConstToReference(parIn, tr, idx)
%TRANSFORMPARCONSTTOREFERENCE Transform nonlinear ROM data.

par = parIn;

T1  = tr.q1.T;
S1  = tr.q1.Tinv;

T2  = tr.q2.T;
S2  = tr.q2.Tinv;

Tx  = tr.qxi.T;
Sx  = tr.qxi.Tinv;

Tg  = tr.qGam.T;

if isfield(par,'Gamma1') && ~isempty(par.Gamma1)
    par.Gamma1 = AeroFlex.sched.transformTensor3(par.Gamma1, T1, S1, S1);
end

if isfield(par,'Gamma2') && ~isempty(par.Gamma2)
    par.Gamma2 = AeroFlex.sched.transformTensor3(par.Gamma2, T1, S2, S2);
end

if isfield(par,'Gamma_g') && ~isempty(par.Gamma_g)
    par.Gamma_g = AeroFlex.sched.transformTensor3(par.Gamma_g, T1, Sx, Sx);
end

if isfield(par,'Gamma_xi') && ~isempty(par.Gamma_xi)
    par.Gamma_xi = AeroFlex.sched.transformTensor3(par.Gamma_xi, Tx, Sx, S1);
end

q1Fields = {'forces_0','N_Thrust','Dw','Ddel','Dddel'};
for k = 1:numel(q1Fields)
    f = q1Fields{k};
    if isfield(par,f) && isnumeric(par.(f)) && size(par.(f),1) == size(T1,2)
        par.(f) = T1*par.(f);
    end
end

gamFields = {'Bw','Bdel','Bddel'};
for k = 1:numel(gamFields)
    f = gamFields{k};
    if isfield(par,f) && isnumeric(par.(f)) && size(par.(f),1) == size(Tg,2)
        par.(f) = Tg*par.(f);
    end
end

par.RateProject = struct('projSet',false,'Pz',[]);
par.compatTransform = true;
end