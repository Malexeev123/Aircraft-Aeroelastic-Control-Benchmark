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
            y = log(round(t_k/cfg.sim.dt)+1) ;
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
    end
end
