classdef PlantRunTime < handle
    properties
        % ---------------- configuration / model handles --------------------
    cfg
    beam
    aero
    base
    idx
    trim

    model                   % AeroFlex.sim.ROMIntegrator

    % ---------------- flexible ROM state -------------------------------
    xFlex                   % current flexible/aeroelastic ROM state
    xFlex0                  % initial flexible/aeroelastic ROM state

    nx double               % flexible ROM state dimension
    nu double               % control input dimension
    nw double               % disturbance dimension

    % ---------------- timing -------------------------------------------
    dt double               % internal plant time step
    t double = 0            % current plant time
    k double = 0            % current plant step index

    % ---------------- body-case / coupled simulation --------------------
    bodyCase string = "wingOnly"
    useRigidBody logical = false

    % ---------------- rigid-body state and parameters -------------------
    rb                      % struct containing rigid-body state
    rbParams                % mass/inertia/geometry/etc.

    % ---------------- flight-condition cache ----------------------------
    flightCond              % struct: Uinf, alpha, beta, euler, omega_B

    % ---------------- last force/moment bookkeeping ---------------------
    last                    % struct storing Clamp6, Ftot_B, Mtot_B, etc.

    % ---------------- trim / fixed rigid inputs -------------------------
    trimThrust double = 0
    trimTailDelta double = 0

    % ---------------- debug --------------------------------------------
    % debugCoupled logical = false
    debugCoupled logical = true
         
    sensor 
    x 
         
         
    
    end
    methods
        function obj = PlantRunTime(cfg, beam, aero, base, x0, trim)
        %PLANTRUNTIME Runtime plant wrapper for wing-only or rigid-flex simulation.
        %
        % Existing call site:
        %   plant = AeroFlex.sim.PlantRunTime(cfg, beam, aero, base, x0)
        %
        % Supports:
        %   cfg.sim.bodyCase = "wing_only"
        %   cfg.sim.bodyCase = "fully_coupled"
        
            obj.cfg  = cfg;
            obj.beam = beam;
            obj.aero = aero;
            obj.base = base;
            obj.trim = trim;
            obj.dt = cfg.sim.dt;
            obj.t  = 0;
            obj.k  = 0;
            obj.x=x0; 
            % --------------------------------------------------------------
            % Body case
            % --------------------------------------------------------------
            if isfield(cfg,'sim') && isfield(cfg.sim,'body_case')
                obj.bodyCase = lower(string(cfg.sim.bodyCase));
            else
                obj.bodyCase = "wingOnly";
            end
        
            obj.useRigidBody = obj.bodyCase == "coupledfull";
            obj.idx = AeroFlex.core.buildIndexStruct(beam.Nm, aero.Na);

            obj.model  = AeroFlex.sim.ROMIntegrator(cfg,beam,aero,base);
            obj.sensor = AeroFlex.sensor.WingVelSensor(beam,cfg);
            obj.nx = size(obj.model.L,1);

            % Failure Handling:
            if isfield(cfg,'nu')
                obj.nu = cfg.nu;
            elseif isfield(cfg,'ctrl') && isfield(cfg.ctrl,'n_surf') && isfield(cfg.ctrl,'var_per')
                obj.nu = cfg.ctrl.n_surf * cfg.ctrl.var_per;
            else
                obj.nu = size(aero.forceMap.B_delta,2) + size(aero.forceMap.B_ddelta,2);
            end
        
            if isfield(cfg,'nw')
                obj.nw = cfg.nw;
            elseif isfield(aero,'forceMap') && isfield(aero.forceMap,'Bw')
                obj.nw = size(aero.forceMap.Bw,2);
            else
                obj.nw = 1;
            end
            % --------------------------------------------------------------
            % Flexible ROM initial state
            % --------------------------------------------------------------
            x0 = x0(:);
        
            if numel(x0) < obj.nx
                error('PlantRunTime:InitialCondition', ...
                      'x0 has length %d, but flexible ROM needs at least nx = %d.', ...
                      numel(x0), obj.nx);
            end
        
            obj.xFlex0 = x0(1:obj.nx);
            obj.xFlex  = obj.xFlex0;
        
            % --------------------------------------------------------------
            % Rigid-body initialization
            % --------------------------------------------------------------
            obj.rbParams = obj.getRigidParams(cfg);
            obj.rb       = obj.initRigidState(cfg, x0);
        
            % --------------------------------------------------------------
            % Fixed trim commands
            % --------------------------------------------------------------
            
            obj.trimThrust = trim.thrust;    

            % if isfield(cfg,'trim') && isfield(cfg.trim,'deltaDeg')
            %     obj.trimTailDelta = deg2rad(cfg.trim.deltaDeg);
            %     if numel(obj.trimTailDelta) > 1
            %         obj.trimTailDelta = obj.trimTailDelta(1);
            %     end
            % else
            %     obj.trimTailDelta = 0;
            % end

            % -------------------------------------------------------------------------
            % Rear elevator/tail trim.
            % Prefer the solved trim struct. Fall back to cfg fields.
            % Stored internally in radians.
            % -------------------------------------------------------------------------
            if nargin >= 6 && isstruct(trim) && isfield(trim,'deltaElev')
                obj.trimTailDelta = trim.deltaElev;
            else
                obj.trimTailDelta = 0;
            end
            
            if numel(obj.trimTailDelta) > 1
                obj.trimTailDelta = obj.trimTailDelta(1);
            end
        
            % --------------------------------------------------------------
            % Flight-condition cache
            % --------------------------------------------------------------
            obj.flightCond = struct();
            obj.flightCond.Uinf    = obj.safeGetFlightSpeed(cfg);
            obj.flightCond.alpha   = obj.safeGetAlpha(cfg);
            obj.flightCond.beta    = 0;
            obj.flightCond.euler   = obj.rb.euler;
            obj.flightCond.omega_B = obj.rb.omega_B;
        
            % --------------------------------------------------------------
            % Last-load cache
            % --------------------------------------------------------------
            obj.last = struct();
            obj.last.Clamp6 = zeros(6,1);
            obj.last.Fwing_B = zeros(3,1);
            obj.last.Mwing_B = zeros(3,1);
            obj.last.Ftail_B = zeros(3,1);
            obj.last.Mtail_B = zeros(3,1);
            obj.last.Ffin_B  = zeros(3,1);
            obj.last.Mfin_B  = zeros(3,1);
            obj.last.Fgrav_B = zeros(3,1);
            obj.last.Fthrust_B = zeros(3,1);
            obj.last.Mthrust_B = zeros(3,1);
            obj.last.Ftot_B = zeros(3,1);
            obj.last.Mtot_B = zeros(3,1);
        
            if isfield(cfg,'debug') && isfield(cfg.debug,'coupledPlant')
                obj.debugCoupled = logical(cfg.debug.coupledPlant);
            else
                obj.debugCoupled = true;
                % obj.debugCoupled = false;
            end
        end
        %--------------------------------------------------------------
        function [z_k,t_k] = readSensors(obj, cfg, log)
            t_k = obj.t;
            % y = log( round(t_k/cfg.sim.dt)+1, :).' ;
            if cfg.forceRealSense
                z_k = obj.sensor.measure(obj.x(obj.model.idx.q1));
            else
                y = log( round(t_k/cfg.sim.dt)+1, :).' ;

                z_k = y;
            end
            % t_k = obj.t;
        end
        %--------------------------------------------------------------
        function x = step(obj,u_k,g_k,mode)
            % [obj.x,~] = obj.model.step(obj.x,u_k,g_k,[],[]);
            switch mode
                case "ROM"          % internal model
                    Nint = round(obj.cfg.ctrl.Ts / obj.cfg.sim.dt);
                    % % if ~exists(S)
                    % 
                    % S = [eye(144), zeros(144, 5)];
                    % % end
                    % x1 = obj.x;
                    for m = 1:Nint
                        [obj.x,~] = obj.model.step(obj.x, u_k, g_k, [], []);
                        % S = [eye(144), zeros(144, 5)];
                        % [obj.x,S] = obj.model.step(obj.x, u_k, g_k,  S, true);
                    end
                    % x2 = obj.x;
                    % % xS = x1+S*[x1; g_k; u_k];
                    % % xS = S*[x1; g_k; u_k];
                    % xS = S(1:length(x1),1:length(x1))*x1;
                    % SensErr = xS-x2;
                    % disp(SensErr);
                    % disp(norm(SensErr));
                    % [obj.x,~] = obj.model.step(obj.x,u_k,g_k,[],[]);
                case "SHARPy"       % read from UDP
                    obj.x = receiveStateFromUDP();              %# your code
                case "openLoop"     % cached log
                    obj.x = getCachedState(obj.t);              %# your code
            end
            % obj.t = obj.t + obj.cfg.sim.dt;
            obj.t = obj.t + obj.cfg.ctrl.Ts;      % ONE sample advance
            x = obj.x;
        end

        function xNext = stepCoupled(obj, u_cmd, gust, rbCmd)
        %STEPCOUPLED One explicit rigid-flexible coupled plant step.
        %
        % Sequence:
        %   1. Read rigid-body flight condition.
        %   2. Update wing flight-condition cache.
        %   3. Propagate flexible wing ROM by one internal dt.
        %   4. Compute wing clamp/root reaction.
        %   5. Compute tail/fin/thrust/gravity loads.
        %   6. Sum forces/moments.
        %   7. Propagate rigid-body equations.
        %   8. Return packed coupled state.
        %
        % This is intentionally explicit coupling. The coupling becomes stiff,
        %% Note to wrap steps 3--7 in a fixed-point subiteration later.
             if nargin < 4 || isempty(rbCmd)
                % rbCmd = struct('delta_e',0,'delta_a',0,'delta_r',0,'thrust',obj.trimThrust);
                rbCmd = [];
             end
            % -------------------------------------------------------------------------
            % Rigid-body command defaults.
            % rbCmd fields are interpreted as increments about trim unless otherwise
            % stated.
            % -------------------------------------------------------------------------
            if isempty(rbCmd)
                rbCmd = struct();
            end
            
            if ~isfield(rbCmd,'delta_e') || isempty(rbCmd.delta_e)
                rbCmd.delta_e = 0;
            end
            
            if ~isfield(rbCmd,'delta_a') || isempty(rbCmd.delta_a)
                rbCmd.delta_a = 0;
            end
            
            if ~isfield(rbCmd,'delta_r') || isempty(rbCmd.delta_r)
                rbCmd.delta_r = 0;
            end
            
            if ~isfield(rbCmd,'thrust') || isempty(rbCmd.thrust)
                rbCmd.thrust = 0;
            end
            
            % Total rear elevator and thrust applied to rigid aircraft.
            delta_e_total = obj.trimTailDelta + rbCmd.delta_e;
            thrust_total  = obj.trimThrust    + rbCmd.thrust;
            rbCmd.delta_e =delta_e_total;
            rbCmd.thrust = thrust_total;
            % --------------------------------------------------------------
            % 0. Validate/sanitize inputs
            % --------------------------------------------------------------
            u_cmd = obj.sanitizeControl(u_cmd);
            gust  = obj.sanitizeDisturbance(gust);
        
            % --------------------------------------------------------------
            % 1. Current rigid-body-derived flight condition
            % --------------------------------------------------------------
            obj.flightCond = obj.computeFlightConditionFromRB();
            
            % --------------------------------------------------------------
            % 2. Push flight condition into flexible ROM, if requested
            % --------------------------------------------------------------
            % gustEff = obj.mapRigidMotionToWingInput(gust);

            obj.updateROMScalingFromRigidBody();
            obj.updateWingFlightCondition();
        
            % --------------------------------------------------------------
            % 3. Propagate flexible wing ROM
            % --------------------------------------------------------------
            obj.xFlex = obj.stepWingROM(u_cmd, gust);
        
            % --------------------------------------------------------------
            % 4. Compute wing clamp/root reaction
            % --------------------------------------------------------------
            Clamp6 = obj.computeWingClampReaction();
            Clamp6 = obj.mirrorWingClamp(Clamp6);

            Fwing_B = Clamp6(1:3);
            Mwing_B = Clamp6(4:6);
        
            % --------------------------------------------------------------
            % 5. Rigid-body aerodynamic / gravity / thrust loads
            % --------------------------------------------------------------
            if isempty(rbCmd )
                [Ftail_B, Mtail_B] = obj.computeTailLoads(u_cmd);
                [Ffin_B,  Mfin_B ] = obj.computeFinLoads(u_cmd);
            
                Fgrav_B = obj.computeGravityBody();
            
                [Fthrust_B, Mthrust_B] = obj.computeThrustLoads(u_cmd);
            else
                % Rigid-body tail/fin/thrust see rbCmd:
                [Ftail_B, Mtail_B] = obj.computeTailLoadsFromCmd(rbCmd);
                [Ffin_B,  Mfin_B ] = obj.computeFinLoadsFromCmd(rbCmd);
                Fgrav_B = obj.computeGravityBody();

                [Fthrust_B, Mthrust_B] = obj.computeThrustLoadsFromCmd(rbCmd);
                
            end
            % --------------------------------------------------------------
            % 6. Sum total aircraft force/moment in body axes
            % --------------------------------------------------------------
            Ftot_B = Fwing_B + Ftail_B + Ffin_B + Fgrav_B + Fthrust_B;
            Mtot_B = Mwing_B + Mtail_B + Mfin_B + Mthrust_B;
        
            % --------------------------------------------------------------
            % 7. Propagate rigid-body EOM
            % --------------------------------------------------------------
            obj.rb = obj.stepRigidBody(obj.rb, Ftot_B, Mtot_B);
        
            % --------------------------------------------------------------
            % 8. Advance time and store diagnostics
            % --------------------------------------------------------------
            obj.k = obj.k + 1;
            obj.t = obj.t + obj.dt;
        
            obj.last.Clamp6 = Clamp6;
            obj.last.Fwing_B = Fwing_B;
            obj.last.Mwing_B = Mwing_B;
            obj.last.Ftail_B = Ftail_B;
            obj.last.Mtail_B = Mtail_B;
            obj.last.Ffin_B  = Ffin_B;
            obj.last.Mfin_B  = Mfin_B;
            obj.last.Fgrav_B = Fgrav_B;
            obj.last.Fthrust_B = Fthrust_B;
            obj.last.Mthrust_B = Mthrust_B;
            obj.last.Ftot_B = Ftot_B;
            obj.last.Mtot_B = Mtot_B;
        
            if obj.debugCoupled
                fprintf(['[stepCoupled] t=%.4f | U=%.3f | alpha=%.3f deg | ', ...
                         '|Fwing|=%.3e | |Ftot|=%.3e | |Mtot|=%.3e\n'], ...
                        obj.t, obj.flightCond.Uinf, rad2deg(obj.flightCond.alpha), ...
                        norm(Fwing_B), norm(Ftot_B), norm(Mtot_B));
            end
        
            % Return packed coupled state.
            xNext = obj.packCoupledState();
        end
        function updateROMScalingFromRigidBody(obj)
        %UPDATEROMSCALINGFROMRIGIDBODY Apply cheap U-dependent scaling.
        %
        % This is only useful if the ROMIntegrator actually uses these fields during
        % model.step. If L already contains the scaling, this does not update L.
        
            mode = "fixed";
            if isfield(obj.cfg,'sim') && isfield(obj.cfg.sim,'rigidToWingMode')
                mode = lower(string(obj.cfg.sim.rigidToWingMode));
            end
        
            if mode ~= "qbar" && mode ~= "alpha_gust"
                return
            end
        
            U0 = obj.cfg.flight.U_inf;
            U  = obj.flightCond.Uinf;
        
            if U0 <= 0
                return
            end
        
            qRatio = (U/U0)^2;
        
            % Store for diagnostics.
            obj.last.qRatio = qRatio;
        
           
        end
    end

    methods (Access = private)
        function xNext = stepWingROM(obj, u_cmd, gust)
        %STEPWINGROM Advance flexible/aeroelastic ROM by one plant time step.
        
            Sdummy = [eye(obj.nx), zeros(obj.nx, obj.nw + obj.nu)];
        
            try
                [xNext, ~] = obj.model.step(obj.xFlex, u_cmd, gust, Sdummy, false, obj.last.qRatio);
            catch
                % Try robustly because weird errors in some commits.
                try
                    [xNext, ~] = obj.model.step(obj.xFlex, u_cmd, gust, Sdummy, true, obj.last.qRatio);
                catch ME
                    error('PlantRunTime:stepWingROM', ...
                          'Could not advance ROMIntegrator.step(). Original error:\n%s', ...
                          ME.message);
                end
            end
        end

        function Clamp6 = computeWingClampReaction(obj)
        %COMPUTEWINGCLAMPREACTION Compute root reaction [Fx;Fy;Fz;Mx;My;Mz].
        %
        % Mirrors trim force path:
        %   xdot = nonlinear_terms(x) + L*x
        %   q1dot = xdot(idx.q1)
        %   Clamp6 = phi1_sA * (Pr*q1dot)
        
            % if ~isfield(obj.beam,'Pr')
            %     error('PlantRunTime:ClampReaction', ...
            %           'beam.Pr is required for clamp reaction extraction.');
            % end
        
            % if ~isfield(obj.beam,'red') || ~isfield(obj.beam.red,'phi1_sA')
            %     error('PlantRunTime:ClampReaction', ...
            %           'beam.red.phi1_sA is required to map clamp modal reactions to 6D resultant.');
            % end
        
            xdot = AeroFlex.sim.nonlinear_terms(obj.xFlex, obj.model.parConst, obj.idx) ...
                 + obj.model.L * obj.xFlex;
        
            q1dot = xdot(obj.idx.q1);
        
            Clamp6 = obj.beam.red.phi1_sA * (obj.beam.Pr * q1dot);
            Clamp6 = obj.mirrorWingClamp(Clamp6);


            Clamp6 = Clamp6(:);
        
            if numel(Clamp6) ~= 6
                error('PlantRunTime:ClampReaction', ...
                      'Clamp6 must have length 6. Got length %d.', numel(Clamp6));
            end
        
            % symmetry treatment if modeling one semi-wing but applying
            % symmetric pair loads to the aircraft.
            if isfield(obj.cfg,'sim') && isfield(obj.cfg.sim,'symmetricWingPair') ...
                    && obj.cfg.sim.symmetricWingPair
        
                % For a symmetric left/right pair:
                %   Fx, Fz, My double.
                %   Fy, Mx, Mz cancel under exact symmetry.
                Clamp6 = [ ...
                    2*Clamp6(1);
                    0;
                    2*Clamp6(3);
                    0;
                    2*Clamp6(5);
                    0 ];
            elseif isfield(obj.cfg,'sim') && isfield(obj.cfg.sim,'wingMultiplicity')
                Clamp6 = obj.cfg.sim.wingMultiplicity * Clamp6;
                % Clamp6 = 2 * Clamp6;
            end
        end

        function fc = computeFlightConditionFromRB(obj)
        %COMPUTEFLIGHTCONDITIONFROMRB Compute Uinf, alpha, beta from body velocity.
        
            vB = obj.rb.v_B(:);
        
            U = norm(vB);
        
            if U < 1e-9
                U = obj.safeGetFlightSpeed(obj.cfg);
                alpha = obj.safeGetAlpha(obj.cfg);
                beta = 0;
            else
                % Body axes convention assumed:
                %   x forward, y right, z down.
                %
                % alpha = atan2(w,u)
                % beta  = asin(v/U)
                alpha = atan2(vB(3), vB(1));
                beta  = asin(max(-1,min(1,vB(2)/U)));
            end
        
            fc = struct();
            fc.Uinf    = U;
            fc.alpha   = alpha;
            fc.beta    = beta;
            fc.euler   = obj.rb.euler(:);
            fc.omega_B = obj.rb.omega_B(:);
        end
        
        function updateWingFlightCondition(obj)
        %UPDATEWINGFLIGHTCONDITION Optionally inject rigid-body attitude/AoA into ROM.
        %
        % By default this is a no-op except for storing obj.flightCond.
        %
        % If cfg.sim.coupleRigidAttitudeIntoWingChi = true and idx.chi exists,
        % then the chi portion of the flexible ROM state is overwritten with
        % a rigid attitude/AoA-compatible value.
        
            if ~isfield(obj.cfg,'sim') || ~isfield(obj.cfg.sim,'coupleRigidAttitudeIntoWingChi')
                return
            end
        
            if ~obj.cfg.sim.coupleRigidAttitudeIntoWingChi
                return
            end
        
            if ~isfield(obj.idx,'chi')
                return
            end
        
            aoaRef = obj.safeGetAlpha(obj.cfg);
        
            % Conservative choice consistent with trim initialization:
            %   chi = [0; alpha - alpha_ref; 0]
            %
            % In works:
            %   chi = obj.flightCond.euler;
            chi = [0;
                   obj.flightCond.alpha - aoaRef;
                   0];
        
            obj.xFlex(obj.idx.chi) = chi;
        end
        function [Ftail_B, Mtail_B] = computeTailLoads(obj, u_cmd)
        %COMPUTETAILLOADS Tail force/moment in body axes.
        
            U     = obj.flightCond.Uinf;
            alpha = obj.flightCond.alpha;
        
            deltaTail = obj.trimTailDelta;
        
            % Later map a surface command to elevator, to do it here.
            if isfield(obj.cfg,'sim') && isfield(obj.cfg.sim,'tailDeltaIndex')
                ii = obj.cfg.sim.tailDeltaIndex;
                if ii >= 1 && ii <= numel(u_cmd)
                    deltaTail = u_cmd(ii);
                end
            end
        
            if exist('tailAeroForceMoment','file') == 2
                [Ftail_B, Mtail_B] = tailAeroForceMoment(U, alpha, deltaTail, obj.cfg);
            else
                Ftail_B = zeros(3,1);
                Mtail_B = zeros(3,1);
            end
        
            Ftail_B = Ftail_B(:);
            Mtail_B = Mtail_B(:);
        end

        function [Ffin_B, Mfin_B] = computeFinLoads(obj, u_cmd)
        %COMPUTEFINLOADS Vertical fin force/moment in body axes.
        
            U    = obj.flightCond.Uinf;
            beta = obj.flightCond.beta;
        
            deltaRudder = 0;
        
            if isfield(obj.cfg,'sim') && isfield(obj.cfg.sim,'rudderIndex')
                ii = obj.cfg.sim.rudderIndex;
                if ii >= 1 && ii <= numel(u_cmd)
                    deltaRudder = u_cmd(ii);
                end
            end
        
            if exist('finAeroForceMoment','file') == 2
                [Ffin_B, Mfin_B] = finAeroForceMoment(U, beta, deltaRudder, obj.cfg);
            else
                Ffin_B = zeros(3,1);
                Mfin_B = zeros(3,1);
            end
        
            Ffin_B = Ffin_B(:);
            Mfin_B = Mfin_B(:);
        end

        function Fgrav_B = computeGravityBody(obj)
        %COMPUTEGRAVITYBODY Gravity force in body axes.
        %
        % Body/NED convention:
        %   inertial gravity vector is [0;0;m*g] in NED.
        %   body axes x forward, y right, z down.
        
            m = obj.rbParams.mass;
            g = obj.rbParams.g;
        
            eul = obj.rb.euler(:);
            C_BI = obj.dcmBodyToInertial(eul);
            C_IB = C_BI.';
        
            Fgrav_I = [0; 0; m*g];
            Fgrav_B = C_IB * Fgrav_I;
        end

        function [Fthrust_B, Mthrust_B] = computeThrustLoads(obj, u_cmd)
        %COMPUTETHRUSTLOADS Thrust force/moment in body axes.
        %
        % Default: constant trim thrust along +x body through CG.
        % TO DO at a later times:
        %   cfg.sim.thrustIndex can map one entry of u_cmd to thrust.
        
            T = obj.trimThrust;
        
            if isfield(obj.cfg,'sim') && isfield(obj.cfg.sim,'thrustIndex')
                ii = obj.cfg.sim.thrustIndex;
                if ii >= 1 && ii <= numel(u_cmd)
                    T = u_cmd(ii);
                end
            end
        
            Fthrust_B = [T; 0; 0];
        
            if isfield(obj.rbParams,'rThrust_B')
                rT = obj.rbParams.rThrust_B(:);
            else
                rT = zeros(3,1);
            end
        
            Mthrust_B = cross(rT, Fthrust_B);
        end

        function rbNext = stepRigidBody(obj, rb, F_B, M_B)
        %STEPRIGIDBODY One explicit RK4 step of rigid-body 6-DoF equations.
        
            dt = obj.dt;
        
            y0 = [rb.r_I(:);
                  rb.v_B(:);
                  rb.euler(:);
                  rb.omega_B(:)];
        
            f = @(y) obj.rigidDerivatives(y, F_B, M_B);
        
            k1 = f(y0);
            k2 = f(y0 + 0.5*dt*k1);
            k3 = f(y0 + 0.5*dt*k2);
            k4 = f(y0 + dt*k3);
        
            y1 = y0 + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
        
            rbNext = rb;
            rbNext.r_I     = y1(1:3);
            rbNext.v_B     = y1(4:6);
            rbNext.euler   = y1(7:9);
            rbNext.omega_B = y1(10:12);
        end

        function ydot = rigidDerivatives(obj, y, F_B, M_B)
        %RIGIDDERIVATIVES Newton-Euler rigid-body equations.
        
            r_I     = y(1:3);
            v_B     = y(4:6);
            euler   = y(7:9);
            omega_B = y(10:12);
        
            r_I_unused = r_I;
        
            m = obj.rbParams.mass;
            I = obj.rbParams.I_B;
        
            C_BI = obj.dcmBodyToInertial(euler);
        
            rdot_I = C_BI * v_B;
        
            % Body-frame translational dynamics:
            %   m*(dv_B/dt + omega x v_B) = F_B
            vdot_B = F_B(:)/m - cross(omega_B, v_B);
        
            % Rotational dynamics:
            %   I*omega_dot + omega x (I*omega) = M_B
            omegaDot_B = I \ (M_B(:) - cross(omega_B, I*omega_B));
        
            eulerDot = obj.eulerRates321(euler, omega_B);
        
            ydot = [rdot_I;
                    vdot_B;
                    eulerDot;
                    omegaDot_B];
        end
        function Clamp6_total = mirrorWingClamp(obj,Clamp6_left)
            S = diag([1,-1,1, -1,1,-1]);
            Clamp6_total = Clamp6_left + S*Clamp6_left;
        end
        function rb = initRigidState(obj, cfg, x0)
        %INITRIGIDSTATE Initialize rigid-body state.
        
            rb = struct();
        
            rb.r_I = zeros(3,1);
        
            U0 = obj.safeGetFlightSpeed(cfg);
            alpha0 = obj.safeGetAlpha(cfg);
        
            % Body-axis velocity. With z down:
            %   u = U*cos(alpha), w = U*sin(alpha)
            rb.v_B = [U0*cos(alpha0);
                      0;
                      U0*sin(alpha0)];
        
            if isfield(cfg,'rigidEOMset') && isfield(cfg.rigidEOMset,'euler0')
                rb.euler = cfg.rigidEOMset.euler0(:);
            elseif isfield(cfg,'trim') && isfield(cfg.trim,'alphaDeg')
                rb.euler = [0; deg2rad(cfg.trim.alphaDeg); 0];
            elseif isfield(cfg,'flight') && isfield(cfg.flight,'aoa_deg')
                rb.euler = [0; deg2rad(cfg.flight.aoa_deg); 0];
            else
                rb.euler = [0; alpha0; 0];
            end
        
            rb.omega_B = zeros(3,1);
        
            % Optional: parse appended rigid states if x0 is longer than flexible nx.
            % Expected appended format:
            %   [r_I(3); v_B(3); euler(3); omega_B(3)]
            if numel(x0) >= obj.nx + 12
                rbTail = x0(obj.nx+1:obj.nx+12);
                rb.r_I     = rbTail(1:3);
                rb.v_B     = rbTail(4:6);
                rb.euler   = rbTail(7:9);
                rb.omega_B = rbTail(10:12);
            end
        end

        function p = getRigidParams(obj, cfg)
        %GETRIGIDPARAMS Retrieve rigid-body mass/inertia parameters.
        
            p = struct();
        
            % ---------------- mass ----------------
            if isfield(cfg,'rigidEOMset') && isfield(cfg.rigidEOMset,'mass')
                p.mass = cfg.rigidEOMset.mass;
            elseif isfield(cfg,'flight') && isfield(cfg.flight,'mass')
                p.mass = cfg.flight.mass;
            else
                % uses this will implement other later.
                try
                    rParams = RigidBody.methods.paramsRigid_PazyUAV();
                    p.mass = rParams.mass;
                catch
                    error('PlantRunTime:RigidParams', ...
                          ['Rigid-body mass not found. Set cfg.rigidEOMset.mass ', ...
                           'or cfg.flight.mass, or make RigidBody.methods.paramsRigid_PazyUAV() available.']);
                end
            end
        
            % ---------------- inertia ----------------
            if isfield(cfg,'rigidEOMset') && isfield(cfg.rigidEOMset,'I_B')
                p.I_B = cfg.rigidEOMset.I_B;
            elseif isfield(cfg,'rigidEOMset') && isfield(cfg.rigidEOMset,'I')
                p.I_B = cfg.rigidEOMset.I;
            else
                try
                    rParams = RigidBody.methods.paramsRigid_PazyUAV();
        
                    if isfield(rParams,'I_B')
                        p.I_B = rParams.I_B;
                    elseif isfield(rParams,'I')
                        p.I_B = rParams.I;
                    elseif all(isfield(rParams,{'Ixx','Iyy','Izz'}))
                        p.I_B = diag([rParams.Ixx, rParams.Iyy, rParams.Izz]);
                    else
                        error('No inertia field found in paramsRigid_PazyUAV.');
                    end
                catch
                    error('PlantRunTime:RigidParams', ...
                          ['Rigid-body inertia not found. Set cfg.rigidEOMset.I_B ', ...
                           'or cfg.rigidEOMset.I.']);
                end
            end
        
            p.I_B = 0.5*(p.I_B + p.I_B.');
        
            if rcond(p.I_B) < 1e-12
                error('PlantRunTime:RigidParams', ...
                      'Rigid-body inertia matrix is singular or ill-conditioned.');
            end
        
            % ---------------- gravity ----------------
            if isfield(cfg,'flight') && isfield(cfg.flight,'g')
                p.g = cfg.flight.g;
            else
                p.g = 9.807;
            end
        
            % ---------------- thrust offset ----------------
            if isfield(cfg,'rigidEOMset') && isfield(cfg.rigidEOMset,'rThrust_B')
                p.rThrust_B = cfg.rigidEOMset.rThrust_B(:);
            else
                p.rThrust_B = zeros(3,1);
            end
        end

        function u = sanitizeControl(obj, u)
        %SANITIZECONTROL Ensure control vector has length nu.
        
            if isempty(u)
                u = zeros(obj.nu,1);
            end
        
            u = u(:);
        
            if isscalar(u) && obj.nu > 1
                u = repmat(u,obj.nu,1);
            end
        
            if numel(u) ~= obj.nu
                error('PlantRunTime:ControlDimension', ...
                      'Control vector has length %d, expected nu = %d.', ...
                      numel(u), obj.nu);
            end
        end

        function w = sanitizeDisturbance(obj, w)
        %SANITIZEDISTURBANCE Ensure disturbance vector has length nw.
        
            if isempty(w)
                w = zeros(obj.nw,1);
            end
        
            w = w(:);
        
            if isscalar(w) && obj.nw > 1
                w = repmat(w,obj.nw,1);
            end
        
            if numel(w) ~= obj.nw
                error('PlantRunTime:DisturbanceDimension', ...
                      'Disturbance vector has length %d, expected nw = %d.', ...
                      numel(w), obj.nw);
            end
        end

        function xPack = packCoupledState(obj)
        %PACKCOUPLEDSTATE Return current plant state.
        %
        % For wing_only:
        %   xPack = xFlex
        %
        % For fully_coupled:
        %   xPack = [xFlex; r_I; v_B; euler; omega_B]
        
            if obj.useRigidBody
                xPack = [obj.xFlex(:);
                         obj.rb.r_I(:);
                         obj.rb.v_B(:);
                         obj.rb.euler(:);
                         obj.rb.omega_B(:)];
            else
                xPack = obj.xFlex(:);
            end
        end

        function C_BI = dcmBodyToInertial(~, euler)
        %DCMBODYTOINERTIAL 3-2-1 body-to-inertial DCM.
        %
        % euler = [phi; theta; psi] = roll, pitch, yaw.
        %
        % Assumes NED-style inertial axes and body axes x-forward, y-right, z-down.
        
            phi   = euler(1);
            theta = euler(2);
            psi   = euler(3);
        
            cphi = cos(phi);   sphi = sin(phi);
            cth  = cos(theta); sth  = sin(theta);
            cpsi = cos(psi);   spsi = sin(psi);
        
            C_BI = [ ...
                cth*cpsi, ...
                sphi*sth*cpsi - cphi*spsi, ...
                cphi*sth*cpsi + sphi*spsi;
        
                cth*spsi, ...
                sphi*sth*spsi + cphi*cpsi, ...
                cphi*sth*spsi - sphi*cpsi;
        
                -sth, ...
                sphi*cth, ...
                cphi*cth ];
        end


        function eulerDot = eulerRates321(~, euler, omega_B)
        %EULERRATES321 3-2-1 Euler kinematics.
        %
        % omega_B = [p; q; r].
        
            phi   = euler(1);
            theta = euler(2);
        
            p = omega_B(1);
            q = omega_B(2);
            r = omega_B(3);
        
            ctheta = cos(theta);
        
            if abs(ctheta) < 1e-8
                warning('PlantRunTime:EulerSingularity', ...
                        'Pitch angle near Euler singularity.');
                ctheta = sign(ctheta)*1e-8;
            end
        
            T = [ ...
                1, sin(phi)*tan(theta),  cos(phi)*tan(theta);
                0, cos(phi),            -sin(phi);
                0, sin(phi)/ctheta,      cos(phi)/ctheta ];
        
            eulerDot = T * [p; q; r];
        end

        function U = safeGetFlightSpeed(~, cfg)
        %SAFEGETFLIGHTSPEED Retrieve nominal freestream speed.
        
            if isfield(cfg,'flight') && isfield(cfg.flight,'U_inf')
                U = cfg.flight.U_inf;
            elseif isfield(cfg,'trim') && isfield(cfg.trim,'U_inf')
                U = cfg.trim.U_inf;
            else
                U = 0;
            end
        end

        function alpha = safeGetAlpha(obj, cfg)
        %SAFEGETALPHA Retrieve nominal/reference angle of attack in radians.
        

            if isprop(obj,'trim') && isfield(obj.trim,'alphaDeg')
                alpha = deg2rad(obj.trim.alphaDeg);
            elseif isfield(cfg,'flight') && isfield(cfg.flight,'aoa_deg')
                alpha = deg2rad(cfg.flight.aoa_deg);
            else
                alpha = 0;
            end
        end
        function [Ftail_B, Mtail_B] = computeTailLoadsFromCmd(obj,rbCmd)
            U = obj.flightCond.Uinf;
            alpha = obj.flightCond.alpha;
        
            if isfield(rbCmd,'delta_e')
                delta_e = rbCmd.delta_e;
            else
                delta_e = 0;
            end
        
            if exist('tailAeroForceMoment','file') == 2
                [Ftail_B,Mtail_B] = tailAeroForceMoment(U,alpha,delta_e,obj.cfg);
            else
                Ftail_B = zeros(3,1);
                Mtail_B = zeros(3,1);
            end
        end
        function [Ffin_B, Mfin_B] = computeFinLoadsFromCmd(obj,rbCmd)
            U = obj.flightCond.Uinf;
            beta = obj.flightCond.beta;
        
            if isfield(rbCmd,'delta_r')
                delta_r = rbCmd.delta_r;
            else
                delta_r = 0;
            end
        
            if exist('finAeroForceMoment','file') == 2
                [Ffin_B,Mfin_B] = finAeroForceMoment(U,beta,delta_r,obj.cfg);
            else
                Ffin_B = zeros(3,1);
                Mfin_B = zeros(3,1);
            end
        end
        function [Fthrust_B, Mthrust_B] = computeThrustLoadsFromCmd(obj,rbCmd)
            if isfield(rbCmd,'thrust')
                T = rbCmd.thrust;
            else
                T = obj.trimThrust;
            end
        
            Fthrust_B = [T;0;0];
        
            if isfield(obj.rbParams,'rThrust_B')
                rT = obj.rbParams.rThrust_B(:);
            else
                rT = zeros(3,1);
            end
        
            Mthrust_B = cross(rT,Fthrust_B);
        end
    end

end
