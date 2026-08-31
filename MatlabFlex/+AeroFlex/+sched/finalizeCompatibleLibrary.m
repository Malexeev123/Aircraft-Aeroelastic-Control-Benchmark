function ROMlib = finalizeCompatibleLibrary(ROMlibIn, varargin)
%FINALIZECOMPATIBLELIBRARY Put ROM points in one reduced coordinate system.

p = inputParser;
p.addParameter('ref_id', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('save_path', '', @(s) ischar(s) || isstring(s));
p.addParameter('allow_identity_aero', false, @islogical);
p.parse(varargin{:});
opt = p.Results;

% ROMlibIn.points.compatibleCoordinates = true;

ROMlib = AeroFlex.sched.loadLibrary(ROMlibIn);

P = ROMlib.points;
Psource = P;
if isempty(opt.ref_id)
    refId = 1;
else
    refId = opt.ref_id;
end

assert(refId >= 1 && refId <= numel(P), 'Invalid reference point id.');

Pref = P(refId);
PrefBasis = Pref.compat.W;

maxErr = struct('q1',0,'qxi',0,'qGam',0);
minRcond = struct('q1',inf,'qxi',inf,'qGam',inf);

for i = 1:numel(P)
    tr = AeroFlex.sched.buildPointTransform(P(i), Pref, ...
        'allow_identity_aero', opt.allow_identity_aero);

    Pnew = AeroFlex.sched.transformPointToReference(P(i), Pref, tr);

    Pnew.compat.W = PrefBasis;
    Pnew.compat.source = Pref.compat.source;
    Pnew.compat.ref_id = refId;
    Pnew.compat.local_to_ref = tr;

    maxErr.q1   = max(maxErr.q1,   tr.err.q1);
    maxErr.qxi  = max(maxErr.qxi,  tr.err.qxi);
    maxErr.qGam = max(maxErr.qGam, tr.err.qGam);

    minRcond.q1   = min(minRcond.q1,   rcond(tr.q1.T));
    minRcond.qxi  = min(minRcond.qxi,  rcond(tr.qxi.T));
    minRcond.qGam = min(minRcond.qGam, rcond(tr.qGam.T));

    Pnew = orderfields(Pnew, P(i));
    P(i) = Pnew;
    % P(i) = Pnew;
end

ROMlib.points = P;
ROMlib.compatibleCoordinates = true;
ROMlib.compatRefId = refId;
ROMlib.compatRefMu = Pref.mu;
ROMlib.compatStats = struct('maxErr',maxErr,'minRcond',minRcond);
ROMlib.method = [char(ROMlib.method) '_compatible'];
ROMlib.validation = AeroFlex.sched.validateCompatibleCoordinates(P,Psource);
ROMlib.physicalRecoveryCompatible = ROMlib.validation.compatible;

fprintf('[finalizeCompatibleLibrary] ref_id=%d mu=[%.6f %.6f]\n', ...
    refId, Pref.mu(1), Pref.mu(2));
fprintf('  max basis residuals: q1=%.3e qxi=%.3e qGam=%.3e\n', ...
    maxErr.q1, maxErr.qxi, maxErr.qGam);
fprintf('  min rcond(T):        q1=%.3e qxi=%.3e qGam=%.3e\n', ...
    minRcond.q1, minRcond.qxi, minRcond.qGam);
fprintf('  transformed checks:  dPz=%.3e dPr=%.3e dPhiVariation=%.3e\n', ...
    ROMlib.validation.dPz, ROMlib.validation.dPr, ROMlib.validation.dPhi);

if ~ROMlib.validation.compatible
    warning('finalizeCompatibleLibrary:CheckFailed', ...
        'Compatible-library conjugacy/recovery validation failed.');
end

if ~isempty(opt.save_path)
    save_path = char(opt.save_path);
    outDir = fileparts(save_path);
    if ~isempty(outDir) && exist(outDir,'dir') ~= 7
        mkdir(outDir);
    end
    save(save_path,'ROMlib','-v7');
end
end
% function ROMlib = finalizeCompatibleLibrary(ROMlibIn, varargin)
% %FINALIZECOMPATIBLELIBRARY Put ROM points in one reduced coordinate system.
%
% p = inputParser;
% p.addParameter('ref_id', [], @(x isempty(x) || isnumeric(x));
% p.addParameter('save_path', '', @(s ischar(s) || isstring(s));
% p.addParameter('allow_identity_aero', false, @islogical);
% p.parse(varargin{:});
% opt = p.Results;
%
% ROMlib = AeroFlex.sched.loadLibrary(ROMlibIn);
% P = ROMlib.points;
%
% if isempty(opt.ref_id)
%     refId = 1;
% else
%     refId = opt.ref_id;
% end
%
% assert(refId >= 1 && refId <= numel(P), 'Invalid reference point id.');
%
% Pref = P(refId);
%
% for i = 1:numel(P)
%     tr = AeroFlex.sched.buildPointTransform(P(i), Pref, ...
%         'allow_identity_aero', opt.allow_identity_aero);
%
%     P(i) = AeroFlex.sched.transformPointToReference(P(i), Pref, tr);
%     P(i).compat = tr;
%     P(i).compat.ref_id = refId;
% end
%
% ROMlib.points = P;
% ROMlib.compatibleCoordinates = true;
% ROMlib.compatRefId = refId;
% ROMlib.compatRefMu = Pref.mu;
% ROMlib.validation = AeroFlex.sched.validateCompatibleCoordinates(P);
%
% if ~ROMlib.validation.compatible
%     warning('finalizeCompatibleLibrary:CheckFailed', ...
%         'Compatible-library validation failed. Check stored transforms.');
% end
%
% if ~isempty(opt.save_path)
%     save_path = char(opt.save_path);
%     outDir = fileparts(save_path);
%     if ~isempty(outDir) && exist(outDir,'dir') ~= 7
%         mkdir(outDir);
%     end
%     save(save_path, 'ROMlib', '-v7');
% end
% end
