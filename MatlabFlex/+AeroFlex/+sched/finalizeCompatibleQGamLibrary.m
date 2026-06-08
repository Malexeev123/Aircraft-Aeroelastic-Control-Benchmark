function ROMlibOut = finalizeCompatibleQGamLibrary(ROMlibIn, savePath)
%FINALIZECOMPATIBLEQGAMLIBRARY Stage-2 coordinate alignment for qGam only.
%
% This version is intentionally conservative: structural coordinates q1,
% q2 and qxi are not transformed, so Gamma tensors remain in their native
% coordinates and are interpolated directly.  The aerodynamic/Krylov qGam
% state can be made compatible if the local aerodynamic basis exists in
% each point.  If no usable basis is available, the input library is returned.

ROMlib = AeroFlex.sched.loadLibrary(ROMlibIn);
ROMlibOut = ROMlib;
ROMlibOut.method = 'stage2_compatible_qGam';
ROMlibOut.compatibility = struct('enabled',false,'reason','No qGam basis found in library points.');

% Try common candidate field names.
basis = cell(numel(ROMlib.points),1);
for i = 1:numel(ROMlib.points)
    p = ROMlib.points(i);
    Zi = [];
    if isfield(p.aero,'Z')
        Zi = p.aero.Z;
    elseif isfield(p.aero,'DataMatrix') && isfield(p.aero.DataMatrix,'Z')
        Zi = p.aero.DataMatrix.Z;
    elseif isfield(p.aero,'ROM') && isfield(p.aero.ROM,'Z')
        Zi = p.aero.ROM.Z;
    end
    if isempty(Zi)
        return
    end
    basis{i} = Zi;
end

Wall = [];
for i = 1:numel(basis), Wall = [Wall, basis{i}]; %#ok<AGROW>
end
[U,~,~] = svd(Wall,'econ');
k = size(basis{1},2);
R = U(:,1:k);

for i = 1:numel(ROMlibOut.points)
    W = basis{i};
    T = R.'*W;
    if rcond(T) < 1e-10
        warning('Point %d qGam transform is ill-conditioned: rcond(T)=%.3e.', i, rcond(T));
    end
    Ti = inv(T);
    idx = ROMlibOut.points(i).idx;
    Nx = ROMlibOut.points(i).dims.Nx;
    Tfull = speye(Nx);
    Tfull(idx.qGam,idx.qGam) = T;
    TfullInv = speye(Nx);
    TfullInv(idx.qGam,idx.qGam) = Ti;

    ROMlibOut.points(i).L = Tfull * ROMlibOut.points(i).L * TfullInv;

    % qGam forcing maps transform with the qGam coordinates.
    pc = ROMlibOut.points(i).parConst;
    pc.Bw    = T * pc.Bw;
    pc.Bdel  = T * pc.Bdel;
    pc.Bddel = T * pc.Bddel;
    ROMlibOut.points(i).parConst = pc;

    ROMlibOut.points(i).compatibility.T_qGam = T;
    ROMlibOut.points(i).compatibility.R_qGam = R;
end

ROMlibOut.compatibility = struct('enabled',true, ...
    'blocks',{{'qGam'}}, ...
    'note','Only qGam was transformed. Structural Gamma tensors are scheduled directly.');

if nargin >= 2 && ~isempty(savePath)
    save(savePath,'ROMlibOut','-v7.3');
    ROMlib = ROMlibOut; %#ok<NASGU>
    save(savePath,'ROMlib','-append');
end
end
