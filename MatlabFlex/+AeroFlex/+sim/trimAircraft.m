function trim = trimAircraft(cfg, beam, aero, base)
% Implements Algorithm 1 (Artola PoC)

% ---------------- constants -------------------------------------------
tolSS  = cfg.trim.tolSteady;      % e.g. 1e‑7
tolNR  = cfg.trim.tolNewton;      % e.g. 1e‑6
% alpha  = cfg.trim.alphaDeg * pi/180;
% alpha  = cfg.trim.alphaDeg * pi/180;

% deltaE = cfg.trim.deltaDeg * pi/180;
deltaE = zeros(4,1) * pi/180;
Tset   = cfg.trim.thrust;         % N
aoa = cfg.flight.aoa_deg* pi/180; 
% disp(aoa)

aoaRef = aoa;
% initial guess will start at linearized val
alpha  = aoaRef;


idx = AeroFlex.core.buildIndexStruct(beam.Nm,aero.Na);
% fThrust = zeros(6,1);
fThrust = 0;
% pre‑compute projectors
Pz  = beam.Pz;                    % 20×20
Pr  = beam.Pr;                    % 20×6

% alpha = alpha - aoa;
% alpha = aoa;
% initial states --------------------------------------------------------
q1  = zeros(beam.Nm,1);
q2  = zeros(beam.Nm,1);
% qxi = zeros(beam.Nm+1,1);
% qxi = zeros(beam.Nm+1,1);  qxi(1) = base.xi_bar(1,1);
qxi = zeros(beam.Nm+1,1); 
qxi(1) = 1;

xi_bar_reord = reshape(base.xi_bar.', [],1);
% disp(norm(xi_bar_reord));
% qxi = base.phi_xi_modes.'*xi_bar_reord

qg  = zeros(cfg.Na,1);
% chi = [0; aoa;0];
chi = [0; alpha-aoaRef; 0];              % Euler [phi, theta, psi] with theta≈aoaRef initially
trim.chi0 = chi;
%%
% % Testing inclusion of the tail Aero effects
% quatRig = eul2quat(chi.');
% xRigid = [quatRig.'; cfg.flight.U_inf;0;0; ... 
%     0;0;0; ...
%     0;0;0];
%%
[phi_xi_trim,   Gamma_xi_iter, Gamma_g_iter, ~, xi_bar_trim] = AeroFlex.core.solve_baseline_quaternion_xi(beam.fem, aero.aeroMesh , ...
beam.phi1, beam.Nm, beam.Mglobal, rad2deg(alpha),beam.phi0);
base.Gamma_xi = Gamma_xi_iter;
base.Gamma_g = Gamma_g_iter;
disp(xi_bar_trim);

x   = [q1;q2;qxi;qg;chi];
sensEq = [eye(size(x)), zeros(size(x,1),1+cfg.ctrl.n_surf*cfg.ctrl.var_per)];
cfg.gust.on      = false;
sim = AeroFlex.sim.SimRunner(cfg,beam,aero,base);

% sim = AeroFlex.sim.ROMIntegrator(cfg,beam,aero,base);
converged = false;  iter = 0;
while ~converged
    iter = iter+1;
    % --- 8‑9 integrate dynamics until Pz*q1 small --------------------
    for k = 1:cfg.trim.itersInt
        [~, x, log_trim, sensEq] = sim.run(x(:,end), cfg.sim.dt, deltaE , cfg.sim.storeSens);   % no gust, no control
        % [x,sensEq] = sim.step(x(:,end),deltaE,0,sensEq,cfg.sim.storeSens);        % no gust
        % disp(size(x));
        q1 = x(sim.idx.q1, end);
        q2 = x(sim.idx.q2, end);
        qxi = x(sim.idx.qxi, end);
        qg = x(sim.idx.qGam, end);
        chi = x(sim.idx.chi, end);
        %% Test
            % qxi = zeros(beam.Nm+1,1); 
            % qxi(1) = 1;
            % x(sim.idx.qxi, end) = qxi;
            % x(sim.idx.chi, end) = zeros(3,1);
        %
        x_nFlex = x(:,end);
        q1 = Pz*q1;   % clamp A‑frame
        x(1:beam.Nm, end) = q1;
        % disp(norm(q1));
        FlexPlant = AeroFlex.sim.ROMIntegrator(cfg,beam,aero,base);
        xNp = AeroFlex.sim.nonlinear_terms(x_nFlex,FlexPlant.parConst,idx);
        xFlex_dot = xNp + FlexPlant.L*x_nFlex;
        % mForces_Sa = beam.Pr*xFlex_dot(idx.q1);
        
        SteadyRes = Pz*xFlex_dot(idx.q1);
        % if cfg.debug.trim
        %     fprintf('#%d  |Steady State Res|=%g\n',iter,norm(SteadyRes));
        % end
        if norm(SteadyRes) < tolSS
            break 
        end
    end
    % log_trim = log_trim(:,end);
    % --- 11 check residual forces at A -------------------------------
    % Res = Pr*q1;               % 3 forces + 3 moments
    Res = Pr*xFlex_dot(1:beam.Nm);               % 3 forces + 3 moments
    q_inf = 0.5*cfg.flight.rho*cfg.flight.U_inf^2;
    a = q_inf*cfg.flight.b_ref^2;

    % Res = a*Res;
    % % % Compute aerodynamic forces and moments from the tail
    % RigidVel = beam.red.phi1_sA*q1;
    % % % Update the residual forces with tail contributions
    % % Res = Res + F_tail;  % Adjust residual forces with tail aerodynamic forces
    % 
    % xRigid(1:4) = eul2quat(chi.').';
    % xRigid(5:7) = RigidVel(1:3);
    % xRigid(8:10) = RigidVel(4:6);
    % [F_tail, M_tail] = tailAeroForceMoment(norm(xRigid(5:7)), alpha, deltaE(1));
    FMwing = beam.red.phi1_sA*Res;
    % FMwing = beam.red.phi1_sA*xFlex_dot(1:beam.Nm);
    % Ftot = FMwing + [F_tail; M_tail];

    if norm(Res) < tolNR
    % if norm(Ftot) < tolNR
        converged = true;
        break
    end
    % Thrust is treated as point force
    thrustLocNodes = 0; % Beam Nodes at which Thrust is performed (need to make this more robust in the config)
    % thrustLocNodes = [2;15]; % Exp of thrust acting on quarter wing half-span
    if isscalar(thrustLocNodes) & thrustLocNodes == 0 % If Thrust is at "clamped" node A
        phi1Thrust = beam.red.phi1_sA;
    else
        phi1Thrust = zeros(6, cfg.Nm );
        for i = 1:length(thrustLocNodes)
            nodeLoc = thrustLocNodes(i);
            nodeDOF = double((nodeLoc-1)*6)+(1:6);
            phi1ThrustLoc = beam.phi1(nodeDOF,:);
            phi1Thrust = phi1Thrust + phi1ThrustLoc;
        end
    end
    JTphi1 = [phi1Thrust.'; zeros(length(x)-cfg.Nm, 6)];
    % uTRig = beam.red.phi1_sA*Tset;
    % xdot = f(cfg, xRigid,uTRig(1:3),Ftot(1:3),Ftot(4:6));
    % --- 14 Newton‑step on α, δe, T ----------------------------------
    % J = beam.dRigdq1;           % 6×20 Jacobian from sensitivities
    J = sensEq(:,:,end);           % 6×20 Jacobian from sensitivities
    
    JwThrust = [J, JTphi1];
    % dU = - (J*Pr)\Res;          % 3×1 update  (α, δe, T)
    % dU = J(idx.q1,:)\(Pr*Res);          % 3×1 update  (α, δe, T)
    dU = J(idx.q1,:)\(Res);          % 3×1 update  (α, δe, T)
    dUwThrust = JwThrust(idx.q1,:)\(Res);
    % dUrig = xdot\xRigid;


    % % Alpha Handling
    % dQuat = base.phiXi_sA*dU(idx.qxi,:);
    % dEul = dU(idx.chi);
    % xi_bar_reg = xi_bar_trim(1,1);
    % dAlpha = quat2eul([xi_bar_reg dQuat(2:end).']); % Recovering q_xi0(t) by setting to 1. Loss of skew-symmetry so check
    % % dAlpha = dEul(2); % Recovering q_xi0(t) by setting to 1. Loss of skew-symmetry so check
    % alpha  = alpha  - dAlpha(2); % Pitch is second
    % 
    % % Delta Handling
    % nu = cfg.ctrl.n_surf*cfg.ctrl.var_per;
    % dDelt = dU(end-nu+1:end); % last 4 is ctrl sens
    % deltaE = deltaE -  dDelt;
%%
    % Alpha Handling
    dQuat = base.phiXi_sA*dUwThrust(idx.qxi,:);
    dEul = dUwThrust(idx.chi);

    xi_bar_reg = xi_bar_trim(1,1);
    dAlpha = quat2eul([xi_bar_reg dQuat(2:end).']); % Recovering q_xi0(t) by setting to 1. Loss of skew-symmetry so check
    % dAlpha = quat2eul(dQuat.'); % Recovering q_xi0(t) by setting to 1. Loss of skew-symmetry so check
    % dAlpha = dEul(2); % Recovering q_xi0(t) by setting to 1. Loss of skew-symmetry so check
    alpha  = alpha  - dAlpha(2); % Pitch is second

    % Delta Handling
    nu = cfg.ctrl.n_surf*cfg.ctrl.var_per;
    gust_start = length(x)+1;
    dDelt = dUwThrust(gust_start+1:gust_start+nu); % 
    deltaE = deltaE -  dDelt;
%%
    % Thrust Handling
    % Tset   = Tset  + dU(3);
    % Tset   = Tset  - dU(idx.q1);
    % fThrust = fThrust - dUwThrust(end-5:end);
    % Tset   = phi1Thrust.'*fThrust; % Temp
    
    deltaThrust = dUwThrust(end-5);
    fThrust = fThrust - deltaThrust; % Fx in Newtons
    Tset   = phi1Thrust.'*[fThrust;0;0;0;0;0]; % 
    % Tset   = Tset; % Temp
    TsetDis = beam.red.phi1_sA*Tset;

    % Update PhiXi info
    [phi_xi_trim,   Gamma_xi_iter, Gamma_g_iter, ~, xi_bar_trim] = AeroFlex.core.solve_baseline_quaternion_xi(beam.fem, aero.aeroMesh , ...
    beam.phi1, beam.Nm, beam.Mglobal, rad2deg(alpha-aoa), beam.phi0);
    % Append Sim
    cfg.trim.thrust = Tset;
    base.Gamma_xi = Gamma_xi_iter;
    base.Gamma_g = Gamma_g_iter;
    sim = AeroFlex.sim.SimRunner(cfg,beam,aero,base);
    % sim = AeroFlex.sim.ROMIntegrator(cfg,beam,aero,base);

    if cfg.debug.trim
        fprintf('#%d  |Trim Force Res|=%g\n',iter,norm(Res));
        % fprintf('#%d  |Trim Force Res|=%g\n',iter,norm(Ftot));
        disp('FMwing'); disp(FMwing);
    end
end

trim.alphaDeg = alpha*180/pi;
trim.deltaDeg = deltaE*180/pi;
trim.thrust   = Tset;
trim.states   = x(:,end);
trim.sensor = log_trim;
end
%%
% % Pull aerodynamic derivatives from LUT instead of single α0
% load(fullfile(AeroFlex.root(),'data','aeroGrid.mat'),'grid','Acell','Bcell');
% [Aero, ~] = AeroFlex.sim.lookupAero(grid,Acell,Bcell,trim.alpha,trim.V); %# << new helper
% trim.Amat = Aero.A;     trim.Bmat = Aero.B;
% disp('trimAircraft: coupled aero matrices inserted from LUT');

function S = skew(v)
S = [   0   -v(3)  v(2);
       v(3)   0   -v(1);
      -v(2)  v(1)   0 ];
end