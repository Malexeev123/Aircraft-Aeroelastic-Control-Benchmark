function trim = TrimRBwFlex(cfg, beam, aero, base)
% TrimRBwFlex — Coupled trim using clamped-wing inner solve (Artola steps 8–10)
% and outer Newton on rigid residuals (Fx, Fz, My) w.r.t. [alpha, delta, T].
%
% Key fixes vs your draft:
%   • PROJECT RATES: enforce q1_dot <- Pz*q1_dot during inner integration.
%   • KEEP RIGID OUT OF INNER LOOP: rigid tail/fin/thrust only in residual.
%   • T as RIGID INPUT: thrust stays scalar in axial residual, not modal forces.
%   • NEWTON JAC: build Jacobian of rigid residual numerically (robust).
%
% Inputs follow your codebase: cfg, beam, aero, base.
% Outputs:
%   trim.states : steady ROM state (x*) at trim
%   trim.alphaDeg, trim.deltaDeg, trim.thrust : trimmed variables
%   trim.sensor : last inner-loop log (if available)
%
% Assumptions:
%   • ROM states/Jacobians are correct (per your note).
%   • beam.Pz (projector), beam.Pr (reaction extractor), beam.red.phi1_sA exist.
%   • cfg.flight has rho, U_inf, mass (or cfg.rigidEOMset.mass).

% ---------------- constants -------------------------------------------
tolSS  = cfg.trim.tolSteady;      % inner steady-state tol on projected q1_dot
tolNR  = cfg.trim.tolNewton;      % outer Newton tol on rigid residual
alpha  = cfg.trim.alphaDeg*pi/180;
% T      = cfg.trim.thrust;         % scalar thrust [N], along +X body
T      = 0;         % scalar thrust [N], along +X body
% elevator (single scalar used by tail model; map your 4-surf vector as needed)
% if isfield(cfg.trim,'deltaDeg'), delta = cfg.trim.deltaDeg*pi/180*ones(4,1); else, delta = 0; end
% if numel(delta) > 1, delta = delta(1); end    % use first entry as elevator
deltaWing = zeros(4,1);
aoaRef = cfg.flight.aoa_deg*pi/180;
% alpha  = alpha - aoaRef;           % Artola: internal states use AoA relative to reference
deltaTa = 0;
deltaElev = deltaTa;

idx  = AeroFlex.core.buildIndexStruct(beam.Nm, aero.Na);
Pz   = beam.Pz;                    % projector on clamped internal motions
Pr   = beam.Pr;                    % extractor of clamp resultants (Nm -> 6)
% 

% % 1a) Pz idempotence
% fprintf('‖Pz^2 - Pz‖ = %.2e\n', norm(Pz*Pz - Pz));
% 
% % 1b) Orthogonality of spaces at A
% rhs = Pr * (Pz * randn(beam.Nm,1));
% fprintf('‖%%A * Pz * v‖ should be ~0: %.2e\n', norm(rhs));
% --- initial ROM state (use your current layout) -----------------------
q1  = zeros(beam.Nm,1);
q2  = zeros(beam.Nm,1);
qxi = zeros(beam.Nm+1,1);  qxi(1) = base.xi_bar(1,1);
qg  = zeros(aero.Na,1);
% chi = [0; aoaRef; 0];              % Euler [phi, theta, psi] with theta≈aoaRef initially
chi = [0; alpha-aoaRef; 0];              % Euler [phi, theta, psi] with theta≈aoaRef initially
x   = [q1; q2; qxi; qg; chi];
trim.chi0 = chi;

% baseline xi coupling for given alpha
[~,Gamma_xi,Gamma_g,~,xi_bar_trim] = ...
    AeroFlex.core.solve_baseline_quaternion_xi(beam.fem, aero.aeroMesh, ...
    beam.phi1, beam.Nm, beam.Mglobal, rad2deg(alpha), beam.phi0);
base.Gamma_xi = Gamma_xi; base.Gamma_g = Gamma_g;

% Build ROM integrator once (we will do explicit qs integration using its L and NL terms)
FlexPlant = AeroFlex.sim.ROMIntegrator(cfg,beam,aero,base);

% mass/weight for residual
% if isfield(cfg,'rigidEOMset') && isfield(cfg.rigidEOMset,'mass')
%     mTOT = cfg.rigidEOMset.mass;
% elseif isfield(cfg.flight,'mass')
%     mTOT = cfg.flight.mass;
% else
%     error('Total mass not provided: set cfg.rigidEOMset.mass or cfg.flight.mass');
% end
rParams = RigidBody.methods.paramsRigid_PazyUAV();
mTOT   = rParams.mass;           % mass of aircraft (kg) - user should set cfg.mass
% mTOT   = rParams.mNoWing/1000;   
g = 9.807; W = mTOT*g;

% outer Newton loop -----------------------------------------------------
maxIt = 100; it = 0;
converged = false;
log_trim = [];
fThrust = 0;
x0 = x;
cfg.trim.useRateProjection = true;
sim = AeroFlex.sim.SimRunner(cfg,beam,aero,base);
% converged = false;  
iter = 0;
while ~converged && it < maxIt
    it = it+1;
    iter = iter+1;
    % x0(idx.chi) = [0;alpha;0];
    % ========== STEPS 8–10: inner clamped solve with RATE PROJECTION ==========
    [x, sim, log_trim] = inner_clamped_steady(x0(:, end), alpha, deltaWing, T, FlexPlant, idx, Pz, tolSS, sim, cfg, beam, aero,base);
    % [x, sim, log_trim] = inner_clamped_steady(x(:, end), alpha, delta, T, FlexPlant, idx, Pz, tolSS, sim, cfg, beam, aero,base);

    % Evaluate clamp forces (A-section resultant) from flexible dynamics
    % Use the *projected* rate at the converged state to form the clamp resultant:
    % xdot = AeroFlex.sim.nonlinear_terms(x(:,end),FlexPlant.parConst,idx) + FlexPlant.L*x(:,end);
    % xdot = AeroFlex.sim.nonlinear_terms(x(:,end),sim.parConst,idx) + FlexPlant.L*x(:,end);
    xdot = AeroFlex.sim.nonlinear_terms(x(:,end),sim.parConst,idx) + sim.L*x(:,end);
    % q1dot_proj = Pz * xdot(idx.q1);
    q1dot_proj = xdot(idx.q1);

    Clamp6 = beam.red.phi1_sA*(Pr*q1dot_proj);   % [Fx;Fy;Fz;Mx;My;Mz] at A
    Clamp6 = mirrorWingClamp(Clamp6);
    % Tail / fin rigid aero (body frame)
    uBody = cfg.flight.U_inf;
    % [F_tail, M_tail] = tailAeroForceMoment(uBody, alpha, delta, cfg);   % user-tunable helper
    [F_tail, M_tail] = tailAeroForceMoment(uBody, alpha, deltaElev, cfg);   % user-tunable helper
    [F_fin,  M_fin ] = finAeroForceMoment (uBody, 0.0  , 0.0 , cfg);    % set to zeros if unused

    % Weight (downwards +Z in NED; in body we take +Z down as well)
    F_w = [0; 0;  W];
    % Thrust (assume through CG along +X)
    F_T = [T; 0; 0]; M_T = [0;0;0];

    % RIGID residual (target zero): choose (Fx, Fz, My)
    F_tot = Clamp6(1:3) + F_tail + F_fin + F_w + F_T;
    M_tot = Clamp6(4:6) + M_tail + M_fin + M_T;
    R = [ F_tot(1); F_tot(3); M_tot(2) ];

    if cfg.debug.trim
        fprintf('It %2d | alpha = %.4f | Wing delta = %.4f | thrust = %.4f | Elevator delta = %.4f  \n',it,alpha,deltaWing(1),T, deltaElev);
        fprintf('|R| = %.3e  [Fx=%.2f, Fz=%.2f, My=%.2f]\n',norm(R),R(1),R(2),R(3));

        % fprintf('It %2d | alpha = %.4f | delta = %.4f | thrust = %.4f \n',it,alpha,delta(1),T);
        % fprintf('|R| = %.3e  [Fx=%.2f, Fz=%.2f, My=%.2f]\n',norm(R),R(1),R(2),R(3));
    end
    if norm(R) < tolNR
        converged = true; break;
    end

    % ================== NEWTON STEP on [alpha, delta, T] ==================
    % Robust finite-difference Jacobian of R w.r.t [alpha, delta, T]
    hA = 1e-3; 
    hE = 1e-3; 
    hD = 1e-3; 
    % hD = 1e-6; 
    hT = max(1e-2, 0.01*abs(T)+1e-3);
    % J = zeros(3,3);
    J = zeros(3,4);
    %%%%%%%%%%%%%%%%%%%%%
    % dR/dalpha
    %%%%%%%%%%%%%%%%%%%%%
    [xA , sim, ~]= inner_clamped_steady(x(:,end), alpha+hA, deltaWing, T,FlexPlant, idx, Pz, tolSS, sim, cfg, beam, aero,base);
    % xdotA = AeroFlex.sim.nonlinear_terms(xA(:,end),FlexPlant.parConst,idx) + FlexPlant.L*xA(:,end);
    % xdotA = AeroFlex.sim.nonlinear_terms(xA(:,end),sim.parConst,idx) + FlexPlant.L*xA(:,end);
    xdotA = AeroFlex.sim.nonlinear_terms(xA(:,end),sim.parConst,idx) + sim.L*xA(:,end);
    % q1dotA = Pz * xdotA(idx.q1);
    q1dotA = xdotA(idx.q1);
    Clamp6A = beam.red.phi1_sA*(Pr*q1dotA);
    Clamp6A = mirrorWingClamp(Clamp6A);

    % [Fta, Mta] = tailAeroForceMoment(uBody, alpha+hA, delta, cfg);
    [Fta, Mta] = tailAeroForceMoment(uBody, alpha+hA, deltaElev, cfg);
    FtotA = Clamp6A(1:3) + Fta + F_fin + F_w + F_T;
    MtotA = Clamp6A(4:6) + Mta + M_fin + M_T;
    JA = [ FtotA(1); FtotA(3); MtotA(2) ];
    J(:,1) = (JA - R)/hA;


    %%%%%%%%%%%%%%%%%%%%%%%%
    % -------------------------------------------------------------------------
    % dR/d deltaElev
    %%%%%%%%%%%%%%%%%%%%
    %  -------------------------------------------------------------------------

    
    [Fte, Mte] = tailAeroForceMoment(uBody, alpha, deltaElev+hE, cfg);
    
    FtotE = Clamp6(1:3) + Fte + F_fin + F_w + F_T;
    MtotE = Clamp6(4:6) + Mte + M_fin + M_T;
    
    RE = [FtotE(1); FtotE(3); MtotE(2)];
    J(:,3) = (RE - R)/hE;
    
    %%%%%%%%%%%%%%%%%%%%%
    % dR/ddelta
    %%%%%%%%%%%%%%%%%%%%%
    [xD , sim, ~] = inner_clamped_steady(x(:,end), alpha, deltaWing+[hD; hD; 0; 0], T, FlexPlant, idx, Pz, tolSS, sim, cfg, beam, aero,base);
    % [xD , sim] = inner_clamped_steady(x(:,end), alpha, delta+[0; 0; hD; hD], T, FlexPlant, idx, Pz, tolSS, sim, cfg, beam, aero,base);
    % xdotD = AeroFlex.sim.nonlinear_terms(xD,FlexPlant.parConst,idx) + FlexPlant.L*xD;
    xdotD = AeroFlex.sim.nonlinear_terms(xD(:,end),sim.parConst,idx) + sim.L*xD(:,end);
    % q1dotD = Pz * xdotD(idx.q1);
    q1dotD = xdotD(idx.q1);
    Clamp6D = beam.red.phi1_sA*(Pr*q1dotD);
    Clamp6D = mirrorWingClamp(Clamp6D);

    % [Ftd, Mtd] = tailAeroForceMoment(uBody, alpha, delta+hD, cfg);
    [Ftd, Mtd] = tailAeroForceMoment(uBody, alpha, deltaElev, cfg);
    % [Ftd, Mtd] = tailAeroForceMoment(uBody, alpha, 0, cfg);
    % [Ftd, Mtd] = tailAeroForceMoment(uBody, alpha, 0, cfg);
    FtotD = Clamp6D(1:3) + Ftd + F_fin + F_w + F_T;
    MtotD = Clamp6D(4:6) + Mtd + M_fin + M_T;
    JD = [ FtotD(1); FtotD(3); MtotD(2) ];
    J(:,2) = (JD - R)/hD;

    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % dR/dT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    JT = [1; 0; 0];  % only Fx residual depends directly on T (no thrust offset -> no pitch moment)
    J(:,4) = JT;
    % dR/dT
    % [xT , sim] = inner_clamped_steady(x(:,end), alpha, delta, T+hT, FlexPlant, idx, Pz, tolSS, sim, cfg, beam, aero,base);
    % xdotT = AeroFlex.sim.nonlinear_terms(xT(:,end),sim.parConst,idx) + sim.L*xT(:,end);
    % % q1dotT = Pz * xdotT(idx.q1);
    % q1dotT = xdotT(idx.q1);
    % Clamp6T = beam.red.phi1_sA*(Pr*q1dotT);
    % Clamp6T = mirrorWingClamp(Clamp6T);
    
    % % [Ftd, Mtd] = tailAeroForceMoment(uBody, alpha, delta+hD, cfg);
    % [FtT, MtT] = tailAeroForceMoment(uBody, alpha, 0, cfg);
    % FtotT = Clamp6T(1:3) + FtT + F_fin + F_w + F_T;
    % MtotT = Clamp6T(4:6) + MtT + M_fin + M_T;
    % JT = [ FtotT(1); FtotT(3); MtotT(2) ];
    % J(:,3) = (JT - R)/hT;
    
    %% Testing Tail Ailerons Trim
    % [Ftda, Mtda] = tailAeroForceMoment(uBody, alpha, deltaTa+hTa, cfg);
    % FtotDa = Clamp6D(1:3) + Ftda + F_fin + F_w + F_T;
    % MtotDa = Clamp6D(4:6) + Mtda + M_fin + M_T;
    % JDa = [ FtotDa(1); FtotDa(3); MtotDa(2) ];
    % % JTmp = (JDa - R)/hTa;
    % J(:,4) = (JDa - R)/hTa;
    


    % Solve Newton system
    % dX = -J\R;

    % if rcond(J) < 1e-10
    %     warning('TrimRBwFlex:IllConditionedNewton', ...
    %         'Trim Newton Jacobian is ill-conditioned. rcond(J)=%.3e. Using pinv.', rcond(J));
    %     dX = -pinv(J)*R;
    % else
        % dX = -J\R;
    % end

    % Least Squares for fourth residual
    Wres = eye(3);
    Jphys = J;
    % Penalize moving wing flap more than elevator/thrust.
    Lambda = diag([ ...
        1e-3, ...   % alpha step regularization
        1e-1,  ...   % wing flap step regularization
        1e-2, ...   % elevator step regularization
        1e-3]);     % thrust step regularization
    
    A = Jphys.'*Wres*Jphys + Lambda;
    b = Jphys.'*Wres*R;
    
    % dX = -A\b;
    if rcond(A) < 1e-12
        warning('TrimRBwFlex:IllConditionedRegLS', ...
            'Regularized trim LS matrix is ill-conditioned. rcond(A)=%.3e. Using pinv.', rcond(A));
        dX = -pinv(A)*b;
    else
        dX = -A\b;
    end

    alpha = alpha + dX(1);
    deltaWing = deltaWing + [dX(2);dX(2);0;0];
    deltaElev = deltaElev + dX(3);
    T     = T     + dX(4);
    % deltaTa = deltaTa + dX(4);
    
    % Update baseline xi coupling after alpha change
    [~,Gamma_xi,Gamma_g,~,xi_bar_trim] = ...
        AeroFlex.core.solve_baseline_quaternion_xi(beam.fem, aero.aeroMesh, ...
        beam.phi1, beam.Nm, beam.Mglobal, rad2deg(alpha),beam.phi0);
    base.Gamma_xi = Gamma_xi; base.Gamma_g = Gamma_g;
    FlexPlant = AeroFlex.sim.ROMIntegrator(cfg,beam,aero,base);  % refresh L, parConst
    
    stateEul = x(idx.chi,end);
    stateEul(2) = alpha-aoaRef;
    % stateEul = [0;alpha-aoaRef;0];
    x(idx.chi,end) = stateEul;
    % x0(idx.chi,end) = stateEul;
    % x(idx.chi,end) = [0;alpha-aoaRef;0];
    x0 = x;
end

% outputs ---------------------------------------------------------------
% trim.alphaDeg = rad2deg(alpha + aoaRef);
trim.alphaDeg = rad2deg(alpha);
% trim.deltaDeg = rad2deg(delta);
trim.deltaDeg = deltaWing;
trim.deltaElev = deltaElev;
trim.thrust   = T;
trim.states   = x(:,end);
trim.sensor   = log_trim;   % (placeholder; keep if you store inner logs)
cfg.trim.useRateProjection = false;
end


function [x, sim, log_trim] = inner_clamped_steady(x0, alpha, delta, T, FlexPlant, idx, Pz, tolSS, sim, cfg, beam, aero,base)
% One quasi-steady pass to converge the clamped wing with q1_dot <- Pz*q1_dot
% Uses explicit Euler on ROM rhs built from the integrator (L and NL terms).

x = x0;
dt = FlexPlant.parConst.dt;
cfg.trim.thrust = beam.red.phi1_sA.'*[T;0;0;0;0;0];
cfg.trim.useRateProjection = true;

[~,Gamma_xi,Gamma_g,~,xi_bar_trim] = ...
    AeroFlex.core.solve_baseline_quaternion_xi(beam.fem, aero.aeroMesh, ...
    beam.phi1, beam.Nm, beam.Mglobal, rad2deg(alpha),beam.phi0);
base.Gamma_xi = Gamma_xi; base.Gamma_g = Gamma_g;
FlexPlant = AeroFlex.sim.ROMIntegrator(cfg,beam,aero,base);

sim = AeroFlex.sim.SimRunner(cfg,beam,aero,base);
% kMax = 10000;
kMax = 20000;
for k = 1:kMax  % cap; your cfg.trim.itersInt can be used instead
    
    [~, x, log_trim, sensEq] = sim.run(x(:,end), cfg.sim.dt, delta, cfg.sim.storeSens);   % no gust, no control
   
            log_trim = log_trim(:,end);
           
            x_nFlex = x(:,end);

            FlexPlant = AeroFlex.sim.ROMIntegrator(cfg,beam,aero,base);
            sim.parConst.u_ctrl = delta;
            xNp = AeroFlex.sim.nonlinear_terms(x_nFlex,sim.parConst,idx);
            xFlex_dot = xNp + sim.L*x_nFlex;

         
            flexRes = Pz*xFlex_dot(idx.q1);
            % flexRes = xFlex_dot(idx.q1);
            Res6 = beam.red.phi1_sA*flexRes;
            %%
            % SteadyRes = Pz*xFlex_dot(idx.q1);
            % SteadyRes = Pz*(xFlex_dot(idx.q1)+ModalR_rb2);
            % SteadyRes = Pz*(xFlex_dot(idx.q1)-ModalR_rb2);
            
            % % if norm(Pz*q1) < tolSS, break, end
            if cfg.debug.trim
            %     % fprintf('#%d  |Coupled Steady State Res|=%g\n',iter,norm(SteadyRes));
                % fprintf('#%d  |Flexible Steady State Res|=%g\n',k,norm(flexRes));
                % fprintf('#%d  |Steady State Res 6x1|=%g\n',k,norm(Res6));
                % disp('|Steady State Res 6x1|');disp(Res6);
            end
            if k == kMax
                % error('Max Steady State Iteration Met')
                warning('Max Steady State Iteration Met')
            end
            % if norm(SteadyRes) < tolSS
            if norm(flexRes) < tolSS
                % fprintf('#%d  |Flexible Steady State Res|=%g\n',k,norm(flexRes));

                break
            end
        
end
end
function [F_tail, M_tail] = tailAeroForceMoment(U, alpha, delta_e, cfg)
% Simple linear tail model around small angles
% if ~isfield(cfg,'tail'), F_tail=[0;0;0]; M_tail=[0;0;0]; return; end
% cfg.tail.S         = 0.45 * cfg.flight.S_ref;   % tail area (assume 20% of wing area)
% cfg.tail.CL_alpha  = 4*pi;                    % lift curve slope (per rad)
%         % cfg.tail.CL_alpha  = 2*pi;                    % lift curve slope (per rad)
% cfg.tail.CL_delta  = 2*pi * 0.8;              % elevator effectiveness (0.8 factor)
%         cfg.tail.CD0       = 0.01;                    % zero-lift drag
%         cfg.tail.e         = 0.8;                     % Oswald efficiency for tail
%         % cfg.tail.x         =  +0.5;                   % tail 0.5 m behind CG
%         cfg.tail.x         =  -0.5;                   % tail 0.5 m behind CG
%         cfg.tail.Z_offset  =  0.0;                    % tail on body axis (no vertical offset)
rho = 1.225;       % air density (kg/m^3) at sea level
S_t = 0.035;       % tail planform area (m^2)
% l_t = -0.6;         % tail moment arm (m, distance aft of CG)
l_t = 0.6;         % tail moment arm (m, distance aft of CG)
a_t = 4.4;         % tail lift curve slope (per radian)
eta = 0.95;        % elevator effectiveness factor (0 < eta <= 1)
cfg.tail.CD0       = 0.01;                    % zero-lift drag
cfg.tail.e         = 0.8;                     % Oswald efficiency for tail
% rho = cfg.flight.rho;
% S_t = cfg.tail.area;   
% CL_a = cfg.tail.CL_alpha;  % per rad
CL_a = a_t;  % per rad
% l_t = cfg.tail.arm;    
% i_t  = cfg.tail.incidence*pi/180;
i_t  = .01*pi/180;
q   = 0.5*rho*U^2;
CLt = CL_a*(alpha - i_t - delta_e);


CD_tail = cfg.tail.CD0 + (CLt^2)/(pi * cfg.tail.e * 4);  % simple induced drag, AR~4 (placeholder)
D_tail = q*S_t * CD_tail;

Lt  = q*S_t*CLt;                 % +Z (down) convention: lift up => negative Z
% F_tail = [0; 0; -Lt];
F_tail = [-D_tail; 0; -Lt];
M_tail = [0; +Lt*l_t; 0];        % nose-up (+My) when tail lift up
end

function [F_fin, M_fin] = finAeroForceMoment(U, beta, rudder, cfg)
% Set to zeros unless you need yaw trim; placeholder included for completeness
F_fin = [0; 0; 0]; M_fin = [0; 0; 0];
end
function Clamp6_total = mirrorWingClamp(Clamp6_left)
    S = diag([1,-1,1, -1,1,-1]);
    Clamp6_total = Clamp6_left + S*Clamp6_left;
end


