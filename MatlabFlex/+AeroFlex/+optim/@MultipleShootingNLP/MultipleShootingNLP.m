classdef MultipleShootingNLP < handle
%MULTIPLESHOOTINGNLP  Builder for direct multiple‑shooting NLP.
%
% Decision vector  p = [x0_2 ... x0_M   u_1 ... u_M   w_1 ... w_M]'
%
% REQUIRED cfg fields (see buildMHEcfg/buildMPCcfg):
%   f,L,Jterm,M,dt,nx,nu,nw,xLow,xUpp,uLow,uUpp,wLow,wUpp
% OPTIONAL:
%   gradFun  – analytic gradient of stage cost
%   jacFun   – analytic Jacobian of dynamics


    properties
        f      % continuous dynamics   xdot = f(x,u,w,t)
        L      % stage cost            L(x,u,w,t)
        Jterm  % terminal cost
        dt     % integration step
        M      % # shooting intervals
        nx, nu, nw
        xLow, xUpp, uLow, uUpp, wLow, wUpp
        stepFcn
        gradFun, jacFun                 % optional exact sensitivities
    end

    methods
        function obj = MultipleShootingNLP(cfg)
            fields = {'f','L','Jterm','dt','M','nx','nu','nw'};
            % fields = fieldnames(cfg);

            for k = 1:numel(fields), obj.(fields{k}) = cfg.(fields{k}); end
            obj.xLow = cfg.xLow; obj.xUpp = cfg.xUpp;
            obj.uLow = cfg.uLow; obj.uUpp = cfg.uUpp;
            obj.wLow = cfg.wLow; obj.wUpp = cfg.wUpp;
            if isfield(cfg,'gradFun'), obj.gradFun = cfg.gradFun; end
            if isfield(cfg,'jacFun' ), obj.jacFun  = cfg.jacFun ; end
        end

        % -------- objective & constraint wrappers for SQP -------------
        function [F,grad_F] = objFun(obj,p,x1)
            % reshape ---------------------------------------------------
            [x0,u,w] = obj.unpack(p);
            xk = x1;           F = 0;   
            if nargout>1, grad_F = zeros(size(p)); end
            % stage costs
            for m = 1:obj.M
                F = F + obj.L(xk,u(:,m),w(:,m));
                if nargout>1 && ~isempty(obj.gradFun)
                    grad_F = grad_F + obj.gradFun(xk,u(:,m),w(:,m),m); 
                end

                xk = xk + obj.f(xk,u(:,m),w(:,m),0)*obj.dt;
                % xk = obj.step(xk,u(:,m),w(:,m));

            end
            % terminal cost

            F = F + obj.Jterm(xk);
        end

        function [c,ceq,Jc,Jceq] = nonlin(obj,p,x1)
            [x0,u,w] = obj.unpack(p);      
            xk = x1;  ceq = [];
            for m = 1:obj.M
                xkNext = xk + obj.f(xk,u(:,m),w(:,m),0)*obj.dt;
                % xkNext = obj.step(xk,u(:,m),w(:,m));

                ceq = [ceq ; x0(:,m) - xk];      
                xk  = xkNext;
            end
            c = [];  Jc = [];  Jceq = [];   % leave Jacobians to SQP fd
        end

        function [x0,u,w] = unpack(obj,p)
            s1 = obj.nx*(obj.M-1);  s2 = s1 + obj.nu*obj.M;
            x0 = reshape(p(1:s1),obj.nx,[]);
            u  = reshape(p(s1+1:s2),obj.nu,[]);
            w  = reshape(p(s2+1:end),obj.nw,[]);
        end

        function p = pack(obj,x0,u,w)
            p = [x0(:) ; u(:) ; w(:)];
        end

        function [lb,ub] = bounds(obj)
            xbox = [obj.xLow obj.xUpp]'; ubox=[obj.uLow obj.uUpp]'; wbox=[obj.wLow obj.wUpp]';
            lb = [ repmat(xbox(1,:),1,obj.M-1) , repmat(ubox(1,:),1,obj.M) , repmat(wbox(1,:),1,obj.M) ]';
            ub = [ repmat(xbox(2,:),1,obj.M-1) , repmat(ubox(2,:),1,obj.M) , repmat(wbox(2,:),1,obj.M) ]';
        end
        function xn = step(obj,x,u,w)
            if isempty(obj.stepFcn)
                xn = x + obj.f(x,u,w,0)*obj.dt;     % default forward Euler
            else
                xn = obj.stepFcn(x,u,w,obj.dt);
            end
        end
    end
end
