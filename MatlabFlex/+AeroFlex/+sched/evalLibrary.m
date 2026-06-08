function sched = evalLibrary(ROMlibIn, mu, cfgLibrary)
%EVALLIBRARY Evaluate/interpolate the scheduled ROM at mu = [U, alpha_deg].
%
% Returned fields are ready to apply to ROMIntegrator/SimRunner/PlantRunTime:
%   sched.L, sched.parConst, sched.idx, sched.beam, sched.base, sched.aero
%
% Stage 1 performs direct componentwise interpolation of all fields consumed
% by nonlinear_terms and the linear IMEX solve.  If a compatible-coordinate
% library has been precomputed, the transformed matrices/fields are used.

if nargin < 3, cfgLibrary = struct(); end
ROMlib = AeroFlex.sched.loadLibrary(ROMlibIn);
[w, ids, info] = AeroFlex.sched.interpWeights(ROMlib, mu, cfgLibrary);
P = ROMlib.points(ids);

sched = struct();
sched.mu = double(mu(:).');
sched.weights = w;
sched.pointIds = ids;
sched.info = info;
sched.method = ROMlib.method;
sched.muNames = ROMlib.muNames;
sched.created = datestr(now,'yyyy-mm-dd HH:MM:SS');

% Linear operator.
vals = cell(numel(P),1);
for k = 1:numel(P), vals{k} = P(k).L; end
sched.L = AeroFlex.sched.lincombNumeric(vals,w);

% Nonlinear bundle: Gamma tensors, force maps, steady loads, scalings.
sched.parConst = AeroFlex.sched.interpParConst(P,w);

% Beam projection/reaction data. These should be identical across the library.
sched.idx = P(1).idx;
sched.beam = P(1).beam;
sched.aero = P(1).aero;
sched.base = P(1).base;

% Baseline fields that vary with alpha.
if isfield(P(1).base,'Gamma_xi')
    vals = arrayfun(@(p){p.base.Gamma_xi},P);
    sched.base.Gamma_xi = AeroFlex.sched.lincombNumeric(vals,w);
end
if isfield(P(1).base,'Gamma_g')
    vals = arrayfun(@(p){p.base.Gamma_g},P);
    sched.base.Gamma_g = AeroFlex.sched.lincombNumeric(vals,w);
end
if isfield(P(1).base,'xi_bar')
    vals = arrayfun(@(p){p.base.xi_bar},P);
    sched.base.xi_bar = AeroFlex.sched.lincombNumeric(vals,w);
end

% Trim/equilibrium offsets.
vals = cell(numel(P),1);
for k = 1:numel(P)
    vals{k} = P(k).x_eq;
end
sched.x_eq = AeroFlex.sched.lincombNumeric(vals,w);

vals = cell(numel(P),1);
for k = 1:numel(P)
    vals{k} = P(k).u_eq;
end
sched.u_eq = AeroFlex.sched.lincombNumeric(vals,w);

sched.trim = P(1).trim;
if isfield(P(1).trim,'alphaDeg')
    vals = arrayfun(@(p){p.trim.alphaDeg},P);
    sched.trim.alphaDeg = AeroFlex.sched.lincombNumeric(vals,w);
end
if isfield(P(1).trim,'deltaDeg')
    vals = arrayfun(@(p){p.trim.deltaDeg},P);
    sched.trim.deltaDeg = AeroFlex.sched.lincombNumeric(vals,w);
end
if isfield(P(1).trim,'deltaElev')
    vals = arrayfun(@(p){p.trim.deltaElev},P);
    sched.trim.deltaElev = AeroFlex.sched.lincombNumeric(vals,w);
end
if isfield(P(1).trim,'thrust')
    vals = arrayfun(@(p){p.trim.thrust},P);
    sched.trim.thrust = AeroFlex.sched.lincombNumeric(vals,w);
end
sched.trim.states = sched.x_eq;

% Diagnostics.
sched.eigL = eig(full(sched.L));
end
