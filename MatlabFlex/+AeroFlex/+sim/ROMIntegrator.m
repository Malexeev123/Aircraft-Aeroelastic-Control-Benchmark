classdef ROMIntegrator < handle
    properties
        L, Ldyn, Lfac, Ufac, piv, gamma, delta, dt, idx, parConst, sched, internalCoupledCoordinate
    end
    methods
        % function obj = ROMIntegrator(cfg,beam,aero,base)
        function obj = ROMIntegrator(cfg,beam,aero,base,sched)
            if nargin<5, sched = []; end
            obj.dt = cfg.sim.dt;
            obj.gamma = (2-sqrt(2))/2;
            obj.delta = -2*sqrt(2)/3;
            obj.idx = AeroFlex.core.buildIndexStruct(beam.Nm,aero.Na);
        
            if nargin >= 5 && ~isempty(sched)
                obj.L = sched.L;
                obj.parConst = sched.parConst;
                obj.idx = sched.idx;
            else
                [obj.L,~] = AeroFlex.core.assemble_L_matrix(cfg,beam,aero,base);
                obj.parConst = AeroFlex.sched.buildParConst(cfg,beam,aero,base);
            end
            
            if ~isfield(obj.parConst,'RateProject') || isempty(obj.parConst.RateProject)
                obj.parConst.RateProject = struct('projSet',false,'Pz',[]);
            end
            if nargin>=5 && ~isempty(sched) && ...
                    isfield(sched,'internalCoupledCoordinate')
                obj.internalCoupledCoordinate=sched.internalCoupledCoordinate;
            else
                obj.internalCoupledCoordinate=struct('enabled',false);
            end
            
            if nargin >= 5 && ~isempty(sched)
                % Scheduled construction must use the same installation
                % contract as the runtime plant. The package owns its step,
                % projection policy, coordinate map, and dependent rebuild.
                obj = AeroFlex.sched.applyToROMIntegrator(obj,sched,cfg);
            else
                obj = obj.rebuildConstrainedOperator();
            end

            % A = speye(size(obj.L)) - obj.gamma*obj.dt*obj.L;
            % [obj.Lfac,obj.Ufac,obj.piv] = lu(A,'vector');
            % obj = obj.rebuildConstrainedOperator();

            % obj.dt = cfg.sim.dt;                           % 1. matrices
            % [obj.L,~] = AeroFlex.core.assemble_L_matrix(cfg,beam,aero,base);
            % obj.gamma = (2-sqrt(2))/2;   obj.delta = -2*sqrt(2)/3;
            % 
            % A = speye(size(obj.L)) - obj.gamma*obj.dt*obj.L;
            % [obj.Lfac,obj.Ufac,obj.piv] = lu(A,'vector');

           
            % obj.parConst = buildParConst(cfg, beam, aero, base);
        end
        %------------------------------------------------------------------
        function [x_np1,S_np1] = step(obj,x_n,u_k,g_k,S_n,storeSens,q_ratio)
        
            if nargin < 6 || isempty(storeSens), storeSens = false; end
            if nargin < 7, q_ratio = []; end
        
            [x_n,S_n]=obj.toInternal(x_n,S_n);
            pc = obj.parConst;
            pc.gust = g_k;
            pc.u_ctrl = u_k;
        
            if ~isempty(q_ratio)
                pc.Fscale    = q_ratio*pc.Fscale;
                pc.scaleA    = q_ratio*pc.scaleA;
                pc.scaleAero = q_ratio*pc.scaleAero;
            end
        
            dt = obj.dt;
            g  = obj.gamma;
            d  = obj.delta;
        
            if isempty(obj.Ldyn)
                obj = obj.rebuildConstrainedOperator();
            end
        
            % Stage 1
            k1N = dt*obj.internalNonlinearTerms(x_n,pc);
            k1N = obj.applyRateProjection(k1N,pc);
        
            if storeSens
                [Nq1,Nu1] = obj.internalNonlinearJacobian(x_n,pc);
                j1N = dt*(Nq1*S_n + Nu1);
                j1N = obj.applyRateProjection(j1N,pc);
        
                rhsS = dt*obj.Ldyn*(S_n + g*j1N);
                j1L = obj.solveImplicit(rhsS);
        
                S1 = S_n + g*(j1L + j1N);
            end
        
            k1L = obj.solveImplicit(dt*obj.Ldyn*(x_n + g*k1N));
            y1  = x_n + g*(k1L + k1N);
        
            % Stage 2
            k2N = dt*obj.internalNonlinearTerms(y1,pc);
            k2N = obj.applyRateProjection(k2N,pc);
        
            if storeSens
                [Nq2,Nu2] = obj.internalNonlinearJacobian(y1,pc);
                j2N = dt*(Nq2*S1 + Nu2);
                j2N = obj.applyRateProjection(j2N,pc);
        
                rhsS = dt*obj.Ldyn*(S_n + (1-g)*j1L + d*j1N + (1-d)*j2N);
                j2L = obj.solveImplicit(rhsS);
        
                S2 = S_n + (1-g)*j1L + g*j2L + d*j1N + (1-d)*j2N;
            end
        
            k2L = obj.solveImplicit(dt*obj.Ldyn*(x_n + (1-g)*k1L + d*k1N + (1-d)*k2N));
        
            % Stage 3
            y2  = x_n + (1-g)*k1L + g*k2L + d*k1N + (1-d)*k2N;
            k3N = dt*obj.internalNonlinearTerms(y2,pc);
            k3N = obj.applyRateProjection(k3N,pc);
        
            if storeSens
                [Nq3,Nu3] = obj.internalNonlinearJacobian(y2,pc);
                j3N = dt*(Nq3*S2 + Nu3);
                j3N = obj.applyRateProjection(j3N,pc);
        
                S_np1 = S_n + (1-g)*(j1L + j2N) + g*(j2L + j3N);
            else
                S_np1 = [];
            end
        
            x_np1 = x_n + (1-g)*(k1L + k2N) + g*(k2L + k3N);
            [x_np1,S_np1]=obj.toPhysical(x_np1,S_np1);
        end
        %------------------------------------------------------------------
        function y = solveImplicit(obj,rhs)
            z = rhs(obj.piv,:);   z = obj.Lfac \ z;   y = obj.Ufac \ z;
        end
        function applySchedule(obj,sched,cfg)
        
            obj.L        = sched.L;
            obj.idx      = sched.idx;
            obj.parConst = sched.parConst;
            if isfield(sched,'internalCoupledCoordinate')
                obj.internalCoupledCoordinate=sched.internalCoupledCoordinate;
            else
                obj.internalCoupledCoordinate=struct('enabled',false);
            end
        
            if isfield(obj.parConst,'dt') && ~isempty(obj.parConst.dt)
                obj.dt = obj.parConst.dt;
            elseif nargin >= 3 && isfield(cfg,'sim') && isfield(cfg.sim,'dt')
                obj.dt = cfg.sim.dt;
            end
        
            if ~isfield(obj.parConst,'RateProject') || isempty(obj.parConst.RateProject)
                obj.parConst.RateProject = struct('projSet',false,'Pz',[]);
            end
        
            obj = obj.rebuildConstrainedOperator();
        end
        function Z = applyRateProjection(obj, Z, pc)
        %APPLYRATEPROJECTION Enforce clamped-root admissible modal rates.
        
            if ~isfield(pc,'RateProject') || ~isstruct(pc.RateProject)
                return
            end
        
            if ~isfield(pc.RateProject,'projSet') || ~pc.RateProject.projSet
                return
            end
        
            if ~isfield(pc.RateProject,'Pz') || isempty(pc.RateProject.Pz)
                return
            end
        
            Pz = pc.RateProject.Pz;
        
            if obj.hasInternalCoupledCoordinate()
                physical=obj.internalCoupledCoordinate.internalToPhysical*Z;
                physical(obj.idx.q1,:)=Pz*physical(obj.idx.q1,:);
                Z=obj.internalCoupledCoordinate.physicalToInternal*physical;
            else
                Z(obj.idx.q1,:) = Pz * Z(obj.idx.q1,:);
            end
        end
        function value = internalNonlinearTerms(obj,state,parConst)
            value=AeroFlex.sched.evaluateInternalCoupledNonlinearTerms( ...
                state,parConst,obj.idx,obj.internalCoupledCoordinate);
        end
        function [Nq,Nu] = internalNonlinearJacobian(obj,state,parConst)
            [Nq,Nu]=AeroFlex.sched.evaluateInternalCoupledNonlinearJacobian( ...
                state,obj.idx,parConst,obj.internalCoupledCoordinate);
        end
        function [state,sensitivity] = toInternal(obj,state,sensitivity)
            if obj.hasInternalCoupledCoordinate()
                state=obj.internalCoupledCoordinate.physicalToInternal*state;
                if ~isempty(sensitivity)
                    sensitivity=obj.internalCoupledCoordinate.physicalToInternal* ...
                        sensitivity;
                end
            end
        end
        function [state,sensitivity] = toPhysical(obj,state,sensitivity)
            if obj.hasInternalCoupledCoordinate()
                state=obj.internalCoupledCoordinate.internalToPhysical*state;
                if ~isempty(sensitivity)
                    sensitivity=obj.internalCoupledCoordinate.internalToPhysical* ...
                        sensitivity;
                end
            end
        end
        function value = hasInternalCoupledCoordinate(obj)
            value=isstruct(obj.internalCoupledCoordinate) && ...
                isfield(obj.internalCoupledCoordinate,'enabled') && ...
                logical(obj.internalCoupledCoordinate.enabled);
        end
        function obj = rebuildConstrainedOperator(obj)
        %REBUILDCONSTRAINEDOPERATOR Rebuild propagated operator and IMEX LU factors.
        %
        % The raw ROM operator obj.L is left untouched.
        %
        % The propagated operator is
        %
        %     Ldyn = Q*L
        %
        % where Q is identity except in the q1 rows:
        %
        %     Ldyn(q1,:) = Pz*L(q1,:)
        %
        % This is the correct operator for the projected continuous-time equation
        %
        %     xdot = Q*(L*x + N(x,u,w)).
        %
        % Consequence:
        %   - project explicit nonlinear increments kN and sensitivity increments jN;
        %   - do NOT project implicit increments kL and jL after solve;
        %   - factor I - gamma*dt*Ldyn, not I - gamma*dt*L.
        
            if isempty(obj.L)
                error('ROMIntegrator:MissingL', ...
                      'Cannot rebuild constrained operator because obj.L is empty.');
            end
        
            if isempty(obj.idx) || ~isfield(obj.idx,'q1')
                error('ROMIntegrator:MissingIndex', ...
                      'Cannot rebuild constrained operator because obj.idx.q1 is missing.');
            end
        
            n = size(obj.L,1);
        
            if size(obj.L,2) ~= n
                error('ROMIntegrator:BadLSize', ...
                      'obj.L must be square. Got %d x %d.', size(obj.L,1), size(obj.L,2));
            end
        
            obj.Ldyn = obj.L;
        
            useProjection = false;
        
            if ~isempty(obj.parConst) && isfield(obj.parConst,'RateProject') && ...
                    isfield(obj.parConst.RateProject,'projSet') && ...
                    logical(obj.parConst.RateProject.projSet)
        
                useProjection = true;
            end
        
            if useProjection
                if ~isfield(obj.parConst.RateProject,'Pz') || isempty(obj.parConst.RateProject.Pz)
                    error('ROMIntegrator:MissingPz', ...
                          'RateProject.projSet is true, but RateProject.Pz is missing.');
                end
        
                Pz = obj.parConst.RateProject.Pz;
                q1 = obj.idx.q1(:);
        
                nq1 = numel(q1);
        
                if ~isequal(size(Pz), [nq1, nq1])
                    error('ROMIntegrator:BadPzSize', ...
                          'Pz must be %d x %d for q1 rows. Got %d x %d.', ...
                          nq1, nq1, size(Pz,1), size(Pz,2));
                end
        
                % Left-projection only. This means Ldyn = Q*L.
                if obj.hasInternalCoupledCoordinate()
                    physicalLdyn=obj.L;
                    physicalLdyn(q1,:)=Pz*physicalLdyn(q1,:);
                    obj.Ldyn=obj.internalCoupledCoordinate.physicalToInternal* ...
                        physicalLdyn*obj.internalCoupledCoordinate.internalToPhysical;
                else
                    obj.Ldyn(q1,:) = Pz*obj.L(q1,:);
                end
            elseif obj.hasInternalCoupledCoordinate()
                obj.Ldyn=obj.internalCoupledCoordinate.Linternal;
            end
        
            A = speye(n) - obj.gamma*obj.dt*obj.Ldyn;
        
            if any(~isfinite(A(:)))
                error('ROMIntegrator:BadImplicitMatrix', ...
                      'Implicit matrix contains NaN or Inf.');
            end
        
            [obj.Lfac, obj.Ufac, obj.piv] = lu(A,'vector');
        end
    end
end
