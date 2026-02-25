function buildAeroGrid(cfg)
% Samples SHARPy at (α, V) grid and stores aerodynamic reduced matrices
%
% OUT: data/aeroGrid.mat
%
alphas = deg2rad(linspace(-2, 10, cfg.Nalpha));   % rad
speeds = linspace(10, 60, cfg.NV);                % m/s
grid   = struct('alpha',alphas,'V',speeds);
Acell  = cell(numel(alphas),numel(speeds));
Bcell  = Acell;
for ia = 1:numel(alphas)
    for iv = 1:numel(speeds)
        [A,B] = runSharpyLin(cfg,alphas(ia),speeds(iv)); %#ok<*NASGU>
        Acell{ia,iv} = A;   Bcell{ia,iv} = B;
        fprintf('lin α=%.1f°  V=%2.0f m/s\n',rad2deg(alphas(ia)),speeds(iv));
    end
end
save(fullfile(cfg.root,'data','aeroGrid.mat'),'grid','Acell','Bcell','-v7.3');
end
