function validateLibrary(ROMlib)
%VALIDATELIBRARY Ensure all points are compatible for interpolation.
assert(isfield(ROMlib,'points') && ~isempty(ROMlib.points), 'ROMlib.points is empty.');
P = ROMlib.points;
ref = P(1).dims;
for k = 1:numel(P)
    assert(P(k).dims.Nm == ref.Nm, 'Point %d has different Nm.', k);
    assert(P(k).dims.Na == ref.Na, 'Point %d has different Na.', k);
    assert(P(k).dims.Nx == ref.Nx, 'Point %d has different state dimension.', k);
    assert(P(k).dims.nu == ref.nu, 'Point %d has different control dimension.', k);
    assert(P(k).dims.nw == ref.nw, 'Point %d has different gust dimension.', k);
    assert(isequal(size(P(k).L), size(P(1).L)), 'Point %d L size mismatch.', k);
    assert(isequal(size(P(k).parConst.Gamma_xi), size(P(1).parConst.Gamma_xi)), 'Point %d Gamma_xi size mismatch.', k);
    assert(isequal(size(P(k).parConst.Gamma_g), size(P(1).parConst.Gamma_g)), 'Point %d Gamma_g size mismatch.', k);
end
end
