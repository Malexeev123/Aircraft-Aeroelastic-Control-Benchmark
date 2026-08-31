function point = extractPoint(cfg, beam, aero, base, trim, varargin)
%EXTRACTPOINT Extract one MATLAB-assembled ROM library point.
%
% The output is the atomic database record used by evalLibrary.  It stores
% both linear and nonlinear quantities because nonlinear_terms depends on
% Gamma1, Gamma2, Gamma_g, Gamma_xi, force maps and steady loads.

p = inputParser;
p.addParameter('body_case','',@(s)ischar(s)||isstring(s));
p.addParameter('source_dir','',@(s)ischar(s)||isstring(s));
p.parse(varargin{:});
opt = p.Results;

idx = AeroFlex.core.buildIndexStruct(beam.Nm,aero.Na);
[L, blk] = AeroFlex.core.assemble_L_matrix(cfg, beam, aero, base);
parConst = AeroFlex.sched.buildParConst(cfg, beam, aero, base);

if nargin < 5 || isempty(trim)
    trim = struct();
end

point = struct();
point.mu = [cfg.flight.U_inf, cfg.flight.aoa_deg];
point.muNames = {'U_inf','alpha_deg'};
% point.name = AeroFlex.sched.pointName(point.mu);
point.name = AeroFlex.sched.pointName(point.mu(1), point.mu(2));
point.created = char(datetime("now","Format","yyyy-MM-dd HH:mm:ss"));
point.source_dir = char(opt.source_dir);
point.body_case = char(opt.body_case);

point.dims = struct();
point.dims.Nm = beam.Nm;
point.dims.Na = aero.Na;
point.dims.Nx = size(L,1);
point.dims.nu = size(aero.forceMap.B_delta,2) + size(aero.forceMap.B_ddelta,2);
point.dims.nw = size(aero.forceMap.Bw,2);

point.cfgFlight = cfg.flight;
point.cfgSim = cfg.sim;
point.cfgCtrl = cfg.ctrl;
point.cfgGust = cfg.gust;

point.L = L;
point.blk = blk;
point.idx = idx;
point.parConst = parConst;

point.base = struct();
point.base.Gamma_xi = base.Gamma_xi;
point.base.Gamma_g  = base.Gamma_g;
point.base.xi_bar   = base.xi_bar;
if isfield(base,'phi_xi_modes'), point.base.phi_xi_modes = base.phi_xi_modes; end
if isfield(base,'FM'), point.base.FM = base.FM; end
if isfield(base,'phiXi_sA'), point.base.phiXi_sA = base.phiXi_sA; end

point.beam = struct();
point.beam.Pz = beam.Pz;
point.beam.Pr = beam.Pr;
point.beam.red = beam.red;
point.beam.eta_e = beam.eta_e;
point.beam.Gamma1 = beam.Gamma1;
point.beam.Gamma2 = beam.Gamma2;
if isfield(beam,'phi0'), point.beam.phi0 = beam.phi0; end
if isfield(beam,'phi1'), point.beam.phi1 = beam.phi1; end
if isfield(beam,'Omega'), point.beam.Omega = beam.Omega; end
if isfield(beam,'Sigma'), point.beam.Sigma = beam.Sigma; end

point.aero = struct();
point.aero.forceMap = aero.forceMap;
point.aero.gust_input = aero.gust_input;
point.aero.forces_aero_beam_dof = aero.DataMatrix.forces_aero_beam_dof;
if isfield(aero,'ROM_dsc'), point.aero.ROM_dsc = aero.ROM_dsc; end
if isfield(aero,'ROM'), point.aero.ROM = aero.ROM; end
if isfield(aero,'DataMatrix'), point.aero.DataMatrix = aero.DataMatrix; end

[Wcompat, Wsrc] = AeroFlex.sched.collectCompatBasis(beam,aero,base);
point.compat = struct();
point.compat.W = Wcompat;
point.compat.source = Wsrc;
point.compat.ref_id = [];
point.compat.local_to_ref = [];

point.trim = trim;
if isfield(trim,'states') && ~isempty(trim.states)
    point.x_eq = trim.states(:);
else
    point.x_eq = zeros(size(L,1),1);
end
if isfield(trim,'u_ctrl') && ~isempty(trim.u_ctrl)
    point.u_eq = trim.u_ctrl(:);
else
    point.u_eq = zeros(point.dims.nu,1);
end

point.validation = struct();
point.validation.eigL = eig(full(L));

point.compatibleCoordinates = true;
end
