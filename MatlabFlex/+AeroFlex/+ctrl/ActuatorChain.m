classdef ActuatorChain < handle
%======================================================================
% ACTUATORCHAIN
%======================================================================
% Implements:
%
%   delta_dot = omega_a * (delta_cmd - delta)
%
% with position and rate saturation.
%
% If var_per == 1:
%   output = delta
%
% If var_per == 2:
%   output = [delta; delta_dot]
%======================================================================

    properties
        Ts double
        omega double

        nSurf double
        varPer double
        nu double

        delta double
        deltaDot double

        deltaL double
        deltaU double
        rateL double
        rateU double
    end

    methods
        function obj = ActuatorChain(cfg, trim)
            obj.Ts = cfg.ctrl.Ts;

            obj.nSurf = cfg.ctrl.n_surf;
            obj.varPer = cfg.ctrl.var_per;
            obj.nu = obj.nSurf*obj.varPer;

            if isfield(cfg,'actuator') && isfield(cfg.actuator,'omega')
                obj.omega = cfg.actuator.omega;
            else
                obj.omega = 40;
            end

            obj.deltaL = obj.expandBound(cfg.uL,obj.nSurf);
            obj.deltaU = obj.expandBound(cfg.uU,obj.nSurf);

            if isfield(cfg,'urateL')
                obj.rateL = obj.expandBound(cfg.urateL,obj.nSurf);
            else
                obj.rateL = -inf(obj.nSurf,1);
            end

            if isfield(cfg,'urateU')
                obj.rateU = obj.expandBound(cfg.urateU,obj.nSurf);
            else
                obj.rateU = inf(obj.nSurf,1);
            end

            if isfield(cfg,'trim') && isfield(cfg.trim,'Uinpt')
                d0 = cfg.trim.Uinpt(:);
            elseif isfield(trim,'Uinpt')
                d0 = trim.Uinpt(:);
            else
                d0 = zeros(obj.nSurf,1);
            end

            if numel(d0) >= obj.nSurf
                d0 = d0(1:obj.nSurf);
            elseif isscalar(d0)
                d0 = repmat(d0,obj.nSurf,1);
            else
                d0 = zeros(obj.nSurf,1);
            end

            obj.delta = min(max(d0,obj.deltaL),obj.deltaU);
            obj.deltaDot = zeros(obj.nSurf,1);
        end

        function [uAct,info] = step(obj,uCmd)
            uCmd = uCmd(:);

            if obj.varPer == 2 && numel(uCmd) == obj.nu
                deltaCmd = uCmd(1:obj.nSurf);
            elseif numel(uCmd) == obj.nSurf
                deltaCmd = uCmd;
            else
                error('ActuatorChain:Dimension', ...
                      'uCmd length %d invalid for nSurf=%d,varPer=%d.', ...
                      numel(uCmd),obj.nSurf,obj.varPer);
            end

            deltaCmdSat = min(max(deltaCmd,obj.deltaL),obj.deltaU);

            deltaDotUnsat = obj.omega * (deltaCmdSat - obj.delta);
            obj.deltaDot = min(max(deltaDotUnsat,obj.rateL),obj.rateU);

            obj.delta = obj.delta + obj.Ts * obj.deltaDot;
            obj.delta = min(max(obj.delta,obj.deltaL),obj.deltaU);

            if obj.varPer == 1
                uAct = obj.delta;
            elseif obj.varPer == 2
                uAct = [obj.delta; obj.deltaDot];
            else
                error('ActuatorChain:Unsupported','Unsupported var_per=%d.',obj.varPer);
            end

            info = struct();
            info.deltaCmd = deltaCmd;
            info.deltaCmdSat = deltaCmdSat;
            info.delta = obj.delta;
            info.deltaDot = obj.deltaDot;
            info.uAct = uAct;
        end

        function [uCmd,info] = previewCommandForEndpoint(obj,uTarget)
        %PREVIEWCOMMANDFORENDPOINT Invert one unsaturated actuator sample.
        %
        % uTarget is the desired next actuator output in the existing
        % [deflection; rate] ordering.  The preview is non-mutating and
        % accepts only targets reproduced by the current first-order actuator
        % without command, position, or rate saturation.
            uTarget = uTarget(:);
            assert(obj.varPer == 2 && numel(uTarget) == obj.nu, ...
                'ActuatorChain:EndpointTarget', ...
                ['Endpoint realization requires a finite [deflection; rate] ', ...
                 'target with the active actuator ordering.']);
            assert(all(isfinite(uTarget)), 'ActuatorChain:EndpointNonfinite', ...
                'Endpoint realization target must be finite.');

            targetDelta = uTarget(1:obj.nSurf);
            targetRate = uTarget(obj.nSurf+1:end);
            commandDelta = obj.delta + targetRate/obj.omega;
            commandRate = targetRate;
            uCmd = [commandDelta; commandRate];

            tolerance = 100 * eps(max(1, norm(uTarget,inf)));
            commandWithinPosition = all(commandDelta >= obj.deltaL-tolerance & ...
                commandDelta <= obj.deltaU+tolerance);
            targetWithinPosition = all(targetDelta >= obj.deltaL-tolerance & ...
                targetDelta <= obj.deltaU+tolerance);
            targetWithinRate = all(targetRate >= obj.rateL-tolerance & ...
                targetRate <= obj.rateU+tolerance);

            commandSat = min(max(commandDelta,obj.deltaL),obj.deltaU);
            rateUnsat = obj.omega*(commandSat-obj.delta);
            realizedRate = min(max(rateUnsat,obj.rateL),obj.rateU);
            realizedDelta = min(max(obj.delta+obj.Ts*realizedRate, ...
                obj.deltaL),obj.deltaU);
            realized = [realizedDelta; realizedRate];
            endpointError = norm(realized-uTarget,inf);
            rateEqualityError = norm(realizedDelta-obj.delta- ...
                obj.Ts*realizedRate,inf);
            accepted = commandWithinPosition && targetWithinPosition && ...
                targetWithinRate && endpointError <= tolerance;

            info = struct('accepted',accepted,'target',uTarget, ...
                'command',uCmd,'realized',realized, ...
                'endpointErrorInfinity',endpointError, ...
                'rateEqualityInfinity',rateEqualityError, ...
                'commandWithinPosition',commandWithinPosition, ...
                'targetWithinPosition',targetWithinPosition, ...
                'targetWithinRate',targetWithinRate, ...
                'tolerance',tolerance);
        end
    end

    methods (Access = private)
        function b = expandBound(~,b,n)
            b = b(:);

            if isscalar(b)
                b = repmat(b,n,1);
            elseif numel(b) ~= n
                error('ActuatorChain:Bound', ...
                      'Bound must be scalar or length %d.',n);
            end
        end
    end
end
