function par = interpParConst(points, w)
%INTERPPARCONST Interpolate the nonlinear parameter bundle.
%
% All fields that enter nonlinear_terms are interpolated.  This is essential:
% Gamma_g and Gamma_xi vary with the baseline quaternion/AoA, while force
% maps and steady loads vary with the SHARPy linearization point.

names = {'Gamma1','Gamma2','Gamma_g','Gamma_xi','forces_0', ...
         'scaleAero','t_inf','Fscale','scaleA','Bw','Dw','Bdel','Ddel', ...
         'Bddel','Dddel','dt','Na','SwTest','N_Thrust'};

par = struct();
for n = 1:numel(names)
    f = names{n};
    if isfield(points(1).parConst,f)
        vals = cell(numel(points),1);
        for k = 1:numel(points)
            vals{k} = points(k).parConst.(f);
        end
        if isnumeric(vals{1}) || islogical(vals{1})
            par.(f) = AeroFlex.sched.lincombNumeric(vals,w);
        else
            par.(f) = vals{1};
        end
    end
end

% Non-numeric bookkeeping should be copied from the closest/nonzero point.
par.gustSet = points(1).parConst.gustSet;
par.gust_input = points(1).parConst.gust_input;
if isfield(points(1).parConst,'RateProject')
    par.RateProject = points(1).parConst.RateProject;
end
end
