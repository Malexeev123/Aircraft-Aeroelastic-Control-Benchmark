function beam = transformBeamToReference(beamIn, beamRef, tr) %#ok<INUSD>
%TRANSFORMBEAMTOREFERENCE Transform stored beam-side ROM fields.

beam = beamIn;

T1 = tr.q1.T;
S1 = tr.q1.Tinv;
S2 = tr.q2.Tinv;

if isfield(beam,'Gamma1') && ~isempty(beam.Gamma1)
    beam.Gamma1 = AeroFlex.sched.transformTensor3(beam.Gamma1,T1,S1,S1);
end

if isfield(beam,'Gamma2') && ~isempty(beam.Gamma2)
    beam.Gamma2 = AeroFlex.sched.transformTensor3(beam.Gamma2,T1,S2,S2);
end

if isfield(beam,'eta_e') && isnumeric(beam.eta_e) && size(beam.eta_e,1) == size(T1,2)
    beam.eta_e = T1*beam.eta_e;
end

if isfield(beam,'Pz'), beam.Pz = T1*beam.Pz*S1; end
if isfield(beam,'Pr'), beam.Pr = T1*beam.Pr*S1; end

if isfield(beam,'red')
    maps = {'phi1_sA','phi2_sA','phi_sA'};
    for k = 1:numel(maps)
        f = maps{k};
        if isfield(beam.red,f) && size(beam.red.(f),2)==size(S1,1)
            beam.red.(f) = beam.red.(f)*S1;
        end
    end
end

if isfield(beam,'phi1'), beam.phi1 = beam.phi1*S1; end
if isfield(beam,'phi0'), beam.phi0 = beam.phi0*S2; end
if isfield(beam,'Omega'), beam.Omega = T1*beam.Omega*S2; end
if isfield(beam,'Sigma'), beam.Sigma = T1*beam.Sigma*S1; end
end
