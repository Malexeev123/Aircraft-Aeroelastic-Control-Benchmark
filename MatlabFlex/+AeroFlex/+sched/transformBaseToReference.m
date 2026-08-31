function base = transformBaseToReference(baseIn, baseRef, tr) %#ok<INUSD>
%TRANSFORMBASETOREFERENCE Transform baseline coupling fields.

base = baseIn;

T1 = tr.q1.T;
S1 = tr.q1.Tinv;

Tx = tr.qxi.T;
Sx = tr.qxi.Tinv;

if isfield(base,'Gamma_g') && ~isempty(base.Gamma_g)
    base.Gamma_g = AeroFlex.sched.transformTensor3(base.Gamma_g,T1,Sx,Sx);
end

if isfield(base,'Gamma_xi') && ~isempty(base.Gamma_xi)
    base.Gamma_xi = AeroFlex.sched.transformTensor3(base.Gamma_xi,Tx,Sx,S1);
end

if isfield(base,'xi_bar') && isnumeric(base.xi_bar) && size(base.xi_bar,1) == size(Tx,2)
    base.xi_bar = Tx*base.xi_bar;
end

if isfield(base,'phi_xi_modes') && size(base.phi_xi_modes,2)==size(Sx,1)
    base.phi_xi_modes = base.phi_xi_modes*Sx;
end

if isfield(base,'phiXi_sA') && size(base.phiXi_sA,2)==size(Sx,1)
    base.phiXi_sA = base.phiXi_sA*Sx;
end
end
