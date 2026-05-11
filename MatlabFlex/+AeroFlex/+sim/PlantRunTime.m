classdef PlantRunTime < handle
    properties, model, sensor, x, t, cfg, end
    methods
        function obj = PlantRunTime(cfg,beam,aero,base,x0)
            obj.model  = AeroFlex.sim.ROMIntegrator(cfg,beam,aero,base);
            obj.sensor = AeroFlex.sensor.WingVelSensor(beam,cfg);
            obj.cfg    = cfg;   obj.x=x0;   obj.t=0;
        end
        %--------------------------------------------------------------
        function [z_k,t_k] = readSensors(obj, cfg, log)
            t_k = obj.t;
            y = log( round(t_k/cfg.sim.dt)+1, :).' ;
            if cfg.forceRealSense
                z_k = obj.sensor.measure(obj.x(obj.model.idx.q1));
            else
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

        function xNext = stepCoupled(obj, u_cmd, gust)
        %STEPCOUPLED One coupled rigid-flexible time step.
        
            % -------------------------------------------------------------
            % 1. Extract current rigid-body state
            % -------------------------------------------------------------
            rb = obj.rbState;
        
            % Compute local flight condition at wing root.
            Uinf  = rb.Uinf;
            alpha = rb.alpha;
            beta  = rb.beta;
            chi   = rb.euler;
            omega = rb.bodyRates;
        
            % -------------------------------------------------------------
            % 2. Update wing ROM parameters from rigid-body motion
            % -------------------------------------------------------------
            obj.updateWingFlightCondition(Uinf, alpha, beta, chi, omega);
        
            % -------------------------------------------------------------
            % 3. Propagate flexible wing ROM
            % -------------------------------------------------------------
            xFlexNext = obj.stepWingROM(u_cmd, gust);
        
            % -------------------------------------------------------------
            % 4. Compute wing reaction/resultant at the fuselage/root
            % -------------------------------------------------------------
            Clamp6 = obj.computeWingClampReaction();
            Fwing_B = Clamp6(1:3);
            Mwing_B = Clamp6(4:6);
        
            % -------------------------------------------------------------
            % 5. Compute rigid-body-only contributions
            % -------------------------------------------------------------
            [Ftail_B, Mtail_B] = obj.computeTailLoads(u_cmd);
            [Ffin_B,  Mfin_B]  = obj.computeFinLoads(u_cmd);
            [Fth_B,   Mth_B]   = obj.computeThrustLoads();
            Fgrav_B            = obj.computeGravityBody();
        
            % -------------------------------------------------------------
            % 6. Sum total aircraft forces/moments
            % -------------------------------------------------------------
            Ftot_B = Fwing_B + Ftail_B + Ffin_B + Fth_B + Fgrav_B;
            Mtot_B = Mwing_B + Mtail_B + Mfin_B + Mth_B;
        
            % -------------------------------------------------------------
            % 7. Propagate rigid-body equations
            % -------------------------------------------------------------
            rbNext = obj.stepRigidBody(Ftot_B, Mtot_B);
        
            % -------------------------------------------------------------
            % 8. Commit states
            % -------------------------------------------------------------
            obj.xFlex  = xFlexNext;
            obj.rbState = rbNext;
        
            xNext = obj.packCoupledState();
        end
    end
end
