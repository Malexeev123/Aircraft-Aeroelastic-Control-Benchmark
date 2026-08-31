function aero = transformAeroToReference(aeroIn, tr)
%TRANSFORMAEROTOREFERENCE Transform stored aerodynamic maps.

aero = aeroIn;

T1 = tr.q1.T;
Tg = tr.qGam.T;

if isfield(aero,'forces_aero_beam_dof') && ...
        isnumeric(aero.forces_aero_beam_dof) && ...
        size(aero.forces_aero_beam_dof,1) == size(T1,2)
    aero.forces_aero_beam_dof = T1*aero.forces_aero_beam_dof;
end

if isfield(aero,'forceMap')
    fm = aero.forceMap;

    q1Fields = {'Dw','D_delta','D_ddelta'};
    for k = 1:numel(q1Fields)
        f = q1Fields{k};
        if isfield(fm,f) && isnumeric(fm.(f)) && size(fm.(f),1) == size(T1,2)
            fm.(f) = T1*fm.(f);
        end
    end

    gamFields = {'Bw','B_delta','B_ddelta'};
    for k = 1:numel(gamFields)
        f = gamFields{k};
        if isfield(fm,f) && isnumeric(fm.(f)) && size(fm.(f),1) == size(Tg,2)
            fm.(f) = Tg*fm.(f);
        end
    end

    aero.forceMap = fm;
end
end