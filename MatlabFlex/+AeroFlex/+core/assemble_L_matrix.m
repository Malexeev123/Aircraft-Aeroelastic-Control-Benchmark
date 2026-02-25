function [L, blk] = assemble_L_matrix(cfg, beam, aero, base)
    %ASSEMBLE_L_MATRIX  Build the coupled linear system for IMEX / trim.
    %
    % Inputs
    %   cfg   – struct from nominalConfig.m       (Nm, Na, etc.)
    %   beam  – BeamModel instance                (phi0, Omega, Dχ …)
    %   aero  – AeroROM instance                  (forceMap with D0, D1, …)
    %   base  – BaselineModes struct              (phi_xi, Gamma_xi, …)
    %
    % Output
    %   L   – sparse or full  (3Nm + 1 + Na + 3) × (...)
    %   blk – struct with block index ranges (R1 … R5)   for debug
    
    Nm = beam.Nm;             Na = cfg.Na;          nChi = 3;
    
    % --- unpack abbreviations ----------------------------------------------
    Omg = beam.Omega;
    I   = eye(Nm);
    Sigma = beam.Sigma;
    FM  = aero.forceMap;
    dt = cfg.sim.dt;
    % disp(base.FM.Bchi)
    % disp(base.FM.Dchi)
    % base.FM.Bchi = eye(Na,3);
    % base.FM.Dchi = eye(Nm,3);
    % R1 = structural momentum + aero coupling ------------------------------

    R1 = [ Sigma + cfg.flight.a*cfg.flight.t_inf*FM.D1 , ...
           Omg - cfg.flight.a*FM.D0/Omg            , ...
           zeros(Nm,Nm+1)                          , ...
           cfg.flight.a*FM.C_Gamma                  , ...
           cfg.flight.a*base.FM.Dchi                    ];
        % R1 = [ Sigma + cfg.flight.a*cfg.flight.t_inf*FM.D1/cfg.flight.Fscale , ...
        %    Omg - cfg.flight.a*FM.D0/cfg.flight.Fscale/Omg            , ...
        %    zeros(Nm,Nm+1)                          , ...
        %    cfg.flight.a*FM.C_Gamma/cfg.flight.Fscale                  , ...
        %    cfg.flight.a*base.FM.Dchi /cfg.flight.Fscale                   ];

      % R1 = [ Sigma + cfg.flight.a*cfg.flight.t_inf^2*FM.D1 , ...
      %      Omg - cfg.flight.a*cfg.flight.t_inf*FM.D0/Omg            , ...
      %      zeros(Nm,Nm+1)                          , ...
      %      cfg.flight.a*cfg.flight.t_inf*FM.C_Gamma                  , ...
      %      cfg.flight.a*cfg.flight.t_inf*base.FM.Dchi                    ];

    % R2 = kinematic block ---------------------------------------------------
    R2 = [ -Omg, zeros(Nm, 2*Nm+1+Na+nChi) ];
    if cfg.struct.discrete
        R2 = dt * R2;                         % scale momentum row
    end
    

    % R3 = strain compatibility (φξ) ----------------------------------------
    R3 = zeros(Nm+1, 3*Nm+1+Na+nChi);
    % R3(:, 2*Nm+1 : 3*Nm+1) = base.phi_xi;     % insert Φξ columns
    
    % R4 = aerodynamic state equation ---------------------------------------
    
    if cfg.aero.discrete
        R4 = [ FM.B1/dt , -FM.B0/Omg/dt , ...
           zeros(Na,Nm+1) , ...
           FM.A_Gamma/dt , base.FM.Bchi/dt ]; 
       
        
        % [ B1_d/dt, -(Ad-I)/dt, zeros(Na,Nm+1), ...
        %        (Ad_gamma-I)/dt,  Bchi_d/dt ];
    else
        % R4 = cfg.flight.t_inf*[ FM.B1 , -FM.B0/Omg/cfg.flight.t_inf , ...
        %    zeros(Na,Nm+1) , ...
        %    FM.A_Gamma/cfg.flight.t_inf , base.FM.Bchi/cfg.flight.t_inf ];
        % 
        R4 = [ FM.B1 , -FM.B0/Omg/cfg.flight.t_inf , ...
           zeros(Na,Nm+1) , ...
           FM.A_Gamma/cfg.flight.t_inf , base.FM.Bchi/cfg.flight.t_inf ];

        % R4 = [ FM.B1 , -FM.B0/Omg, ...
        %    zeros(Na,Nm+1) , ...
        %    FM.A_Gamma/cfg.flight.t_inf , base.FM.Bchi/cfg.flight.t_inf ];
    end
    % R5 = rigid‑body χ̂ eqs (simple rate form) ------------------------------
    % R5 = [ FM.X, zeros(nChi, 2*Nm+1+Na+nChi) ]; % From Prev EXP.
    R5 = [ base.FM.X, zeros(nChi, 2*Nm+1+Na+nChi) ]; % Calculated
    

    if cfg.struct.retain_rigid
        R1 = R1 + [ Kss , Csr*0 , zeros(Nf,Nf+1) , zeros(Nf,Na) , Csr ];
        R5 = R5 + [ Krs , zeros(6,Nf) , zeros(6,Nf+1+Na) , Crr ];
    end

    % --- assemble -----------------------------------------------------------
    L = [ R1 ; R2 ; R3 ; R4 ; R5 ];
    
    
    % optional block index struct
    blk.R1 = 1:Nm;
    blk.R2 = Nm+1 : 2*Nm;
    blk.R3 = 2*Nm+1 : 3*Nm+1;
    blk.R4 = 3*Nm+2 : 3*Nm+1+Na;
    blk.R5 = 3*Nm+1+Na+1 : size(L,1);

    %% Validation
    % % build continuous & discrete L
    % [Lc,~] = assemble_L_matrix(cfg,beam,aero,base);           % aero.continuous
    % cfg.aero.discrete = true;
    % [Ld,~] = assemble_L_matrix(cfg,beam,aero,base);           % aero.discrete
    % 
    % % Rayleigh quotient → modal damping ζi ≈ –Re(λ)/|λ|
    % lam_c = eig(Lc);   damp_c = -real(lam_c)./abs(lam_c);
    % lam_d = eig(Ld);   damp_d = -log(abs(lam_d))/dt ./ abs(angle(lam_d)/dt);
    % 
    % assert(all(abs(damp_c-damp_d) < 1e-3), ...
    %        'Damping mismatch >1e‑3 for Δt/t_inf ≤ 0.1')
    % 
    % 
    % tags = {'force','momentum','strain','circulation','rate'};
    % blk  = {'R1','R2','R3','R4','R5'};
    % for k = 1:5
    %     row = Lc.(blk{k});
    %     test = sqrt(sum(row.^2,2));          % each eqn magnitude
    %     assert(all(isfinite(test)),'Units broke in %s',blk{k})
    % end
    % disp('Unit sanity passed')

    % 
    % cfg.aero.discrete = false;  [Lc,~] = assemble_L_matrix(cfg,beam,aero,base);
    % cfg.aero.discrete = true;   [Ld,~] = assemble_L_matrix(cfg,beam,aero,base);
    % 
    % Δt   = cfg.time.dt;         % 6.25e‑4
    % lamC = eig(Lc);
    % dampC= -real(lamC)./abs(lamC);
    % 
    % lamD = eig(Ld);
    % wnD  = angle(lamD)/Δt;          % natural freq
    % dampD= -log(abs(lamD))./sqrt(log(abs(lamD)).^2 + wnD.^2);
    % 
    % assert( max(abs(dampC-dampD)) < 1e‑3, 'scaling mismatch' )

end

%% No sub terms
% Define system matrices for linear part (Eq 30)
% L_matrix = [Sigma,   Omega, zeros(num_modes,num_modes+1), aero_scale*C_Gamma, aero_scale*D_chi;
%            -Omega,                zeros(num_modes),           zeros(num_modes,num_modes+1+num_aero_states2+3);
%            zeros(num_modes+1, 3*num_modes+1+num_aero_states2+3);
%            B1,            (-B0/Omega),                 zeros(num_aero_states2,num_modes+1), A_Gamma_scaled, B_chi/t_inf;
%            X,                     zeros(3,num_modes),        zeros(3,num_modes+1+num_aero_states2+3)];

% % Full
% if scale_Lmat == 1
% % Define system matrices for linear part (Eq 30)
% L_matrix = [Sigma + aero_scale*t_inf*D1,   Omega-aero_scale*D0/Omega, zeros(num_modes,num_modes+1), aero_scale*C_Gamma, aero_scale*D_chi;
%            -Omega,                zeros(num_modes),           zeros(num_modes,num_modes+1+num_aero_states2+3);
%            zeros(num_modes+1, 3*num_modes+1+num_aero_states2+3);
%            B1,            (-B0/Omega),                 zeros(num_aero_states2,num_modes+1), A_Gamma_scaled, B_chi/t_inf;
%            X,                     zeros(3,num_modes),        zeros(3,num_modes+1+num_aero_states2+3)];
% 
% else
%     % No scaling
%     % Define system matrices for linear part (Eq 30)
%     L_matrix = [Sigma + D1,   Omega-D0/Omega, zeros(num_modes,num_modes+1), C_Gamma, D_chi;
%                -Omega,                zeros(num_modes),           zeros(num_modes,num_modes+1+num_aero_states2+3);
%                zeros(num_modes+1, 3*num_modes+1+num_aero_states2+3);
%                B1,            (-B0/Omega),                 zeros(num_aero_states2,num_modes+1), A_Gamma, B_chi;
%                X,                     zeros(3,num_modes),        zeros(3,num_modes+1+num_aero_states2+3)];
% 
%      % 
%      % L_matrix = [Sigma + aero_scale*D1,   Omega-aero_scale*D0/Omega, zeros(num_modes,num_modes+1), aero_scale*C_Gamma, aero_scale*D_chi;
%      %           -Omega,                zeros(num_modes),           zeros(num_modes,num_modes+1+num_aero_states2+3);
%      %           zeros(num_modes+1, 3*num_modes+1+num_aero_states2+3);
%      %           B1,            (-B0/Omega),                 zeros(num_aero_states2,num_modes+1), A_Gamma, B_chi;
%      %           X,                     zeros(3,num_modes),        zeros(3,num_modes+1+num_aero_states2+3)];
% end