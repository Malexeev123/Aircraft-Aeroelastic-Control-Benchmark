function [W,src] = collectCompatBasis(beam,aero,base)
%COLLECTCOMPATBASIS Collect bases for coordinate alignment.

Nm  = beam.Nm;
Na  = aero.Na;
Nxi = Nm + 1;

[Wq1,  s1] = getFirstExisting(beam, { ...
    {'phi1'}, ...
    {'red','phi1'}, ...
    {'red','Phi1'}, ...
    {'red','ModeVars_continuous','phi1'}, ...
    {'red','ModeVars_continuous','Phi1'}, ...
    {'red','ModeVars_discrete','phi1'}, ...
    {'red','ModeVars_discrete','Phi1'} });

[Wxi, sxi] = getFirstExisting(base, { ...
    {'phi_xi_modes'}, ...
    {'Phi_xi_modes'}, ...
    {'phi_xi'}, ...
    {'Phi_xi'} });
aero.DataMatrix.KrylovBasis = aero.DataMatrix.aero_gamma_state.project_V.r.r;
[Wg, sg] = getFirstExisting(aero, { ...
    {'DataMatrix','KrylovBasis'}, ...
    {'DataMatrix','Z'}, ...
    {'DataMatrix','V'}, ...
    {'DataMatrix','W'}, ...
    {'DataMatrix','KrylovSharpyContSS','UserData','Z'}, ...
    {'DataMatrix','KrylovSharpyContSS','UserData','V'}, ...
    {'DataMatrix','KrylovSharpyContSS','UserData','W'} });
    % {'ROM_dsc','UserData','Z'}, ...
    % {'ROM_dsc','UserData','V'}, ...
    % {'ROM_dsc','UserData','W'}, ...
    % {'ROM_cts','UserData','Z'}, ...
    % {'ROM_cts','UserData','V'}, ...
    % {'ROM_cts','UserData','W'} });
        % {'DataMatrix','aero_gamma_state', 'project_V'}, ...


if isempty(Wg)
    [Wg, sg] = findAeroBasisInH5(aero, Na);
end

if isempty(Wq1)
    error('collectCompatBasis:MissingQ1', ...
        'No q1 basis found. Store beam.phi1 or beam.red.ModeVars_continuous.phi1.');
end

if isempty(Wxi)
    error('collectCompatBasis:MissingQxi', ...
        'No qxi basis found. Store base.phi_xi_modes.');
end

if isempty(Wg)
    error('collectCompatBasis:MissingQGam', ...
        ['No qGam/Krylov basis found. Store it in aero.DataMatrix.KrylovBasis, ', ...
         'ROM_dsc.UserData.V/Z/W, or the HDF5 projection file.']);
end

W.q1  = orientBasis(Wq1, Nm,  'q1');
W.qxi = orientBasis(Wxi, Nxi, 'qxi');

% SHARPy can store a leading gust-monitoring block in the Krylov basis.
% That block is not part of the MATLAB reduced qGam state and its size can
% vary with flight condition.  Strip it only for basis matching.
WgRaw = orientBasis(Wg, Na, 'qGam');
[W.qGam, qGamStrip] = stripQGamGustRows(WgRaw, aero);

src.q1   = s1;
src.qxi  = sxi;
src.qGam = sg;

% Keep these as diagnostics.  They are useful when checking the raw library.
src.qGamRawRows       = qGamStrip.rawRows;
src.qGamRows          = qGamStrip.finalRows;
src.qGamStrippedRows  = qGamStrip.nStrip;
src.qGamVsizeID       = qGamStrip.VsizeID;
src.qGamStripApplied  = qGamStrip.stripped;
end


% -------------------------------------------------------------------------
function [val,src] = getFirstExisting(S,paths)

val = [];
src = '';

for i = 1:numel(paths)
    p = paths{i};
    tmp = S;
    ok = true;

    for j = 1:numel(p)
        name = p{j};

        if ~hasMember(tmp,name)
            ok = false;
            break
        end

        tmp = getMember(tmp,name);

        if isempty(tmp)
            ok = false;
            break
        end
    end

    if ok
        val = tmp;
        src = strjoin(p,'.');
        return
    end
end
end

% -------------------------------------------------------------------------
function tf = hasMember(S,name)

tf = false;

if isstruct(S)
    tf = isfield(S,name);
elseif isobject(S)
    % if ~isempty(S.(name))
        tf = isprop(S,name);
    % end
end
end

% -------------------------------------------------------------------------
function val = getMember(S,name)

if isstruct(S)
    val = S.(name);
elseif isobject(S) && isprop(S,name)
    val = S.(name);
else
    val = [];
end
end

% -------------------------------------------------------------------------
function W = orientBasis(W,ncol,name)

W = full(double(W));

if ~ismatrix(W)
    error('collectCompatBasis:BadBasisDims', ...
        '%s basis must be 2-D. Got ndims=%d.', name, ndims(W));
end

if size(W,2) == ncol
    return
end

if size(W,1) == ncol
    W = W.';
    return
end

error('collectCompatBasis:BadBasisSize', ...
    '%s basis must have %d columns. Got %d x %d.', ...
    name, ncol, size(W,1), size(W,2));
end

% -------------------------------------------------------------------------
function [W,src] = findAeroBasisInH5(aero,Na)

W = [];
src = '';

h5Fields = {'vproj_h5','rom_h5','data_h5'};

for k = 1:numel(h5Fields)
    f = h5Fields{k};

    if ~hasMember(aero,f)
        continue
    end

    h5file = getMember(aero,f);

    if isempty(h5file)
        continue
    end

    h5file = char(h5file);

    if exist(h5file,'file') ~= 2
        continue
    end

    try
        info = h5info(h5file);
        [Wcand, scand] = scanH5Group(h5file, info, Na);
    catch
        Wcand = [];
        scand = '';
    end

    if ~isempty(Wcand)
        W = Wcand;
        src = sprintf('%s:%s', f, scand);
        return
    end
end
end

% -------------------------------------------------------------------------
function [bestW,bestSrc] = scanH5Group(h5file,grp,Na)

bestW = [];
bestSrc = '';
bestScore = -inf;

% Datasets in this group.
for i = 1:numel(grp.Datasets)
    ds = grp.Datasets(i);
    dpath = joinH5Path(grp.Name, ds.Name);

    sz = ds.Dataspace.Size;
    if numel(sz) ~= 2
        continue
    end

    if ~(any(sz == Na) && max(sz) > Na)
        continue
    end

    lname = lower(dpath);

    score = 0;
    if contains(lname,'krylov'), score = score + 5; end
    if contains(lname,'basis'),  score = score + 5; end
    if contains(lname,'proj'),   score = score + 4; end
    if contains(lname,'vproj'),  score = score + 4; end
    if contains(lname,'right'),  score = score + 2; end
    if contains(lname,'v'),      score = score + 1; end
    if contains(lname,'w'),      score = score + 1; end

    if score <= bestScore
        continue
    end

    try
        A = h5read(h5file,dpath);
    catch
        continue
    end

    if ~isnumeric(A) || isempty(A) || ~isreal(A)
        continue
    end

    if ~(any(size(A) == Na) && max(size(A)) > Na)
        continue
    end

    bestW = A;
    bestSrc = dpath;
    bestScore = score;
end

% Recurse into subgroups.
for i = 1:numel(grp.Groups)
    [Wsub, Ssub] = scanH5Group(h5file, grp.Groups(i), Na);

    if isempty(Wsub)
        continue
    end

    lsub = lower(Ssub);
    score = 0;
    if contains(lsub,'krylov'), score = score + 5; end
    if contains(lsub,'basis'),  score = score + 5; end
    if contains(lsub,'proj'),   score = score + 4; end
    if contains(lsub,'vproj'),  score = score + 4; end
    if contains(lsub,'right'),  score = score + 2; end
    if contains(lsub,'v'),      score = score + 1; end
    if contains(lsub,'w'),      score = score + 1; end

    if score > bestScore
        bestW = Wsub;
        bestSrc = Ssub;
        bestScore = score;
    end
end
end

% -------------------------------------------------------------------------
function p = joinH5Path(g,name)

if strcmp(g,'/')
    p = ['/' name];
else
    p = [g '/' name];
end
end

% -------------------------------------------------------------------------
function [Wout,info] = stripQGamGustRows(Win,aero)
%STRIPQGAMGUSTROWS Remove SHARPy's variable gust block from the qGam basis.
%
% This is used only for Stage-2 coordinate compatibility.
% It must not be applied to ROM_cts, forceMap, parConst, Bw, Dw, or states.
%
% Expected DataMatrix.aero_gamma_state.VsizeID:
%   no gust-state block : [64; 640; 64; 64]
%   gust-state included : [ng; 64; 640; 64; 64]
%
% The first block, ng, is a gust-monitoring bookkeeping block.  Its size can
% vary with flight condition.  The remaining rows are the common aerodynamic
% parent-state layout used for qGam basis comparison.

Wout = Win;

info = struct();
info.rawRows  = size(Win,1);
info.finalRows = size(Win,1);
info.nStrip = 0;
info.stripped = false;
info.VsizeID = [];

if ~hasMember(aero,'DataMatrix')
    return
end

DM = getMember(aero,'DataMatrix');

if ~hasMember(DM,'aero_gamma_state')
    return
end

ags = getMember(DM,'aero_gamma_state');

if ~hasMember(ags,'VsizeID')
    return
end

VsizeID = double(getMember(ags,'VsizeID'));
VsizeID = VsizeID(:);

info.VsizeID = VsizeID;

if numel(VsizeID) == 4
    expectedRows = sum(VsizeID);

    if size(Win,1) ~= expectedRows
        error(['collectCompatBasis:BadQGamVsizeID\n', ...
               'qGam basis has %d rows, but sum(VsizeID)=%d.'], ...
               size(Win,1), expectedRows);
    end

    return
end

if numel(VsizeID) ~= 5
    error(['collectCompatBasis:UnknownQGamVsizeID\n', ...
           'Expected VsizeID length 4 or 5. Got length %d.'], ...
           numel(VsizeID));
end

nStrip = VsizeID(1);
expectedRawRows = sum(VsizeID);
expectedFinalRows = sum(VsizeID(2:end));

if size(Win,1) == expectedFinalRows
    % Already stripped upstream.  Leave it alone.
    info.finalRows = size(Win,1);
    return
end

if size(Win,1) ~= expectedRawRows
    error(['collectCompatBasis:BadQGamVsizeID\n', ...
           'qGam basis row count does not match VsizeID.\n', ...
           'size(W,1)=%d, sum(VsizeID)=%d, sum(VsizeID(2:end))=%d.'], ...
           size(Win,1), expectedRawRows, expectedFinalRows);
end

Wout = Win(nStrip+1:end,:);

info.nStrip = nStrip;
info.stripped = true;
info.finalRows = size(Wout,1);

if info.finalRows ~= expectedFinalRows
    error(['collectCompatBasis:BadQGamStrip\n', ...
           'Unexpected qGam row count after stripping. Got %d, expected %d.'], ...
           info.finalRows, expectedFinalRows);
end
end
% function [W,src] = collectCompatBasis(beam,aero,base)
% %COLLECTCOMPATBASIS Collect reduced bases for library alignment.
% 
% Nm  = beam.Nm;
% Na  = aero.Na;
% Nxi = Nm + 1;
% 
% [Wq1,  s1] = getFirstExisting(beam, { ...
%     {'phi1'}, ...
%     {'red','phi1'}, ...
%     {'red','Phi1'}, ...
%     {'red','ModeVars_continuous','phi1'}, ...
%     {'red','ModeVars_continuous','Phi1'}, ...
%     {'red','ModeVars_discrete','phi1'}, ...
%     {'red','ModeVars_discrete','Phi1'} });
% 
% [Wxi, sxi] = getFirstExisting(base, { ...
%     {'phi_xi_modes'}, ...
%     {'Phi_xi_modes'}, ...
%     {'phi_xi'}, ...
%     {'Phi_xi'} });
% 
% [Wg, sg] = getFirstExisting(aero, { ...
%     {'compat','W','qGam'}, ...
%     {'ROM','Z'}, ...
%     {'ROM','V'}, ...
%     {'ROM','basis','Z'}, ...
%     {'ROM_dsc','Z'}, ...
%     {'ROM_dsc','V'}, ...
%     {'ROM_dsc','basis','Z'}, ...
%     {'DataMatrix','Z'}, ...
%     {'DataMatrix','V'}, ...
%     {'DataMatrix','KrylovBasis'}, ...
%     {'DataMatrix','basis','Z'} });
% 
% if isempty(Wq1)
%     error('collectCompatBasis:MissingQ1', ...
%         'No q1 basis found. Store beam.phi1 or beam.red.ModeVars_continuous.phi1.');
% end
% 
% if isempty(Wxi)
%     error('collectCompatBasis:MissingQxi', ...
%         'No qxi basis found. Store base.phi_xi_modes.');
% end
% 
% if isempty(Wg)
%     error('collectCompatBasis:MissingQGam', ...
%         'No qGam/Krylov basis found. Store aero.ROM.Z, aero.ROM_dsc.Z, or aero.DataMatrix.Z.');
% end
% 
% W.q1   = orientBasis(Wq1, Nm,  'q1');
% W.qxi  = orientBasis(Wxi, Nxi, 'qxi');
% W.qGam = orientBasis(Wg,  Na,  'qGam');
% 
% src.q1   = s1;
% src.qxi  = sxi;
% src.qGam = sg;
% end
% 
% function W = orientBasis(W,ncol,name)
% W = full(W);
% 
% if size(W,2) == ncol
%     return
% end
% 
% if size(W,1) == ncol
%     W = W.';
%     return
% end
% 
% error('collectCompatBasis:BadBasisSize', ...
%     '%s basis must have %d columns. Got %d x %d.', ...
%     name, ncol, size(W,1), size(W,2));
% end
% 
% function [val,src] = getFirstExisting(S,paths)
% val = [];
% src = '';
% 
% for i = 1:numel(paths)
%     p = paths{i};
%     tmp = S;
%     ok = true;
% 
%     for j = 1:numel(p)
%         if ~(isstruct(tmp)||isa(tmp, 'AeroFlex.beam.BeamModel') ||isa(tmp, 'AeroFlex.aero.AeroROM')) || ...
%             ~(isfield(tmp,p{j})|| isnumeric(tmp.(p{j}))) || isempty(tmp.(p{j}))
%             ok = false;
%             break
%         end
%         tmp = tmp.(p{j});
%     end
% 
%     if ok
%         val = tmp;
%         src = strjoin(p,'.');
%         return
%     end
% end
% end
