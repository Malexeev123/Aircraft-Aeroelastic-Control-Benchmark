classdef ComplementaryCommandFusion < handle
%======================================================================
% COMPLEMENTARYCOMMANDFUSION
%======================================================================
% Fuses slow outer-loop command with fast inner-loop nMPC command:
%
%   u = uTrim + LP(uOuter-uTrim) + HP(uInner-uTrim)
%
% where:
%   LP(s) = wc/(s+wc)
%   HP(s) = s/(s+wc) = 1 - LP(s)
%
% Discrete implementation:
%   yLP[k] = a*yLP[k-1] + (1-a)*x[k],
%   a = exp(-wc*Ts)
%======================================================================

    properties
        Ts double
        wc double
        a double

        uTrim double

        yOuterLP double
        yInnerLP double

        initialized logical = false
        innerIsDeviation logical = false
    end

    methods
        function obj = ComplementaryCommandFusion(cfg, trim)
            obj.Ts = cfg.ctrl.Ts;

            if isfield(cfg,'fusion') && isfield(cfg.fusion,'omega_c')
                obj.wc = cfg.fusion.omega_c;
            else
                obj.wc = 2.5;
            end

            obj.a = exp(-obj.wc*obj.Ts);

            nSurf = cfg.ctrl.n_surf;
            varPer = cfg.ctrl.var_per;
            nu = nSurf*varPer;

            if isfield(cfg,'trim') && isfield(cfg.trim,'Uinpt')
                u0 = cfg.trim.Uinpt(:);
            elseif isfield(trim,'uTrim')
                u0 = trim.uTrim(:);
            else
                u0 = zeros(nu,1);
            end

            if numel(u0) == nSurf && varPer == 2
                u0 = [u0; zeros(nSurf,1)];
            elseif isscalar(u0)
                u0 = repmat(u0,nu,1);
            elseif numel(u0) ~= nu
                u0 = zeros(nu,1);
            end

            obj.uTrim = u0(:);

            obj.yOuterLP = zeros(nu,1);
            obj.yInnerLP = zeros(nu,1);

            if isfield(cfg,'fusion') && isfield(cfg.fusion,'innerIsDeviation')
                obj.innerIsDeviation = logical(cfg.fusion.innerIsDeviation);
            end
        end

        function [uCmd,info] = fuse(obj,uOuter,uInner)
            uOuter = uOuter(:);
            uInner = uInner(:);

            if numel(uOuter) ~= numel(obj.uTrim)
                error('ComplementaryCommandFusion:Dimension', ...
                      'uOuter length %d does not match uTrim length %d.', ...
                      numel(uOuter), numel(obj.uTrim));
            end

            if numel(uInner) ~= numel(obj.uTrim)
                error('ComplementaryCommandFusion:Dimension', ...
                      'uInner length %d does not match uTrim length %d.', ...
                      numel(uInner), numel(obj.uTrim));
            end

            dOuter = uOuter - obj.uTrim;

            if obj.innerIsDeviation
                dInner = uInner;
            else
                dInner = uInner - obj.uTrim;
            end

            obj.yOuterLP = obj.a*obj.yOuterLP + (1-obj.a)*dOuter;
            obj.yInnerLP = obj.a*obj.yInnerLP + (1-obj.a)*dInner;

            dOuterLow = obj.yOuterLP;
            dInnerHigh = dInner - obj.yInnerLP;

            uCmd = obj.uTrim + dOuterLow + dInnerHigh;

            info = struct();
            info.uOuter = uOuter;
            info.uInner = uInner;
            info.dOuterLow = dOuterLow;
            info.dInnerHigh = dInnerHigh;
            info.uCmd = uCmd;
            info.wc = obj.wc;
        end

        function uBase = previewOuterBaseHorizon(obj,uOuter,nHorizon)
        %PREVIEWOUTERBASEHORIZON Non-mutating low-pass outer-command preview.
            assert(isscalar(nHorizon) && isfinite(nHorizon) && ...
                nHorizon >= 1 && nHorizon == round(nHorizon), ...
                'ComplementaryCommandFusion:Horizon', ...
                'nHorizon must be a positive integer.');
            if isvector(uOuter)
                uOuter = repmat(obj.validateCommand(uOuter,'uOuter'),1,nHorizon);
            else
                assert(isnumeric(uOuter) && isreal(uOuter) && ...
                    isequal(size(uOuter),[numel(obj.uTrim),nHorizon]) && ...
                    all(isfinite(uOuter),'all'), ...
                    'ComplementaryCommandFusion:OuterPreviewSequence', ...
                    ['A varying outer-command preview must have one finite ', ...
                     'command column per horizon interval.']);
            end

            yOuter = obj.yOuterLP;
            uBase = zeros(numel(obj.uTrim),nHorizon);
            for j = 1:nHorizon
                dOuter = uOuter(:,j) - obj.uTrim;
                yOuter = obj.a*yOuter + (1-obj.a)*dOuter;
                uBase(:,j) = obj.uTrim + yOuter;
            end
        end

        function [uInner,info] = rawInnerForApplied(obj,uApplied,uOuter)
        %RAWINNERFORAPPLIED Invert the next fusion step without mutation.
            uApplied = obj.validateCommand(uApplied,'uApplied');
            uOuter = obj.validateCommand(uOuter,'uOuter');
            assert(isfinite(obj.a) && obj.a > 0 && obj.a <= 1, ...
                'ComplementaryCommandFusion:Coefficient', ...
                'The fusion coefficient must lie in (0,1].');

            uBase = obj.previewOuterBaseHorizon(uOuter,1);
            dInner = obj.yInnerLP + (uApplied-uBase(:,1))/obj.a;
            if obj.innerIsDeviation
                uInner = dInner;
            else
                uInner = obj.uTrim + dInner;
            end
            assert(all(isfinite(uInner)), ...
                'ComplementaryCommandFusion:NonfiniteInverse', ...
                'The reconstructed raw inner command must be finite.');

            info = struct('uApplied',uApplied,'uOuter',uOuter, ...
                'uBase',uBase(:,1),'uInner',uInner);
        end
    end

    methods(Access=private)
        function u = validateCommand(obj,u,name)
            u = u(:);
            assert(numel(u) == numel(obj.uTrim), ...
                'ComplementaryCommandFusion:Dimension', ...
                '%s length %d does not match uTrim length %d.', ...
                name,numel(u),numel(obj.uTrim));
            assert(all(isfinite(u)), 'ComplementaryCommandFusion:Nonfinite', ...
                '%s must contain only finite values.',name);
        end
    end
end
