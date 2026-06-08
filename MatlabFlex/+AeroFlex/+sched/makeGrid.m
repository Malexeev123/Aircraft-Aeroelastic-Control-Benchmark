function grid = makeGrid(cfg)
%MAKEGRID Build the U-alpha library grid from cfg.library.
cfg = AeroFlex.sched.defaultLibraryConfig(cfg);
U = cfg.library.U_grid(:).';
a = cfg.library.alpha_grid_deg(:).';
[UU,AA] = ndgrid(U,a);
mu = [UU(:), AA(:)];

names = cell(size(mu,1),1);
for k = 1:size(mu,1)
    names{k} = AeroFlex.sched.pointName(mu(k,:));
end

grid = struct();
grid.mu = mu;
grid.U_grid = U;
grid.alpha_grid_deg = a;
grid.names = names;
grid.muNames = cfg.library.muNames;
end
