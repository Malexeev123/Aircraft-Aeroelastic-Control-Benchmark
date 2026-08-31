function Pout = transformPointToReference(P, Pref, tr)
%TRANSFORMPOINTTOREFERENCE Transform one ROM point to reference coordinates.

Pout = P;

Pout.L = tr.T * P.L * tr.Tinv;

if isfield(P,'x_eq') && ~isempty(P.x_eq)
    Pout.x_eq = tr.T * P.x_eq(:);
end

if isfield(P,'trim') && isfield(P.trim,'states') && ~isempty(P.trim.states)
    Pout.trim.states = tr.T * P.trim.states(:);
end

Pout.parConst = AeroFlex.sched.transformParConstToReference(P.parConst, tr, P.idx);

if isfield(P,'beam')
    Pout.beam = AeroFlex.sched.transformBeamToReference(P.beam, Pref.beam, tr);
end

if isfield(P,'aero')
    Pout.aero = AeroFlex.sched.transformAeroToReference(P.aero, tr);
end

if isfield(P,'base')
    Pout.base = AeroFlex.sched.transformBaseToReference(P.base, Pref.base, tr);
end

Pout.idx = Pref.idx;
Pout.compatibleCoordinates = true;
end