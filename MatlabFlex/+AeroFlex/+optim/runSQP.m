function [p,f,H,info] = runSQP(nlp,p,H,opts)
%RUNSQP  Sequential Quadratic Programming with BFGS & bound projection
%
%   [p,H,info] = RUNSQP(nlp,p,H,opts)
%
%   Inputs
%   ------
%   nlp.cost     : @(p) [f,g]                 (objective and gradient)
%   nlp.nonl     : @(p) [ceq,Aeq]             (eq rows and Jacobian)
%   nlp.lb, ub   : box bounds  (vectors)
%
%   p            : warm-start   (column)
%   H            : SPD Hessian approximation (same size as p)
%   opts fields  : MaxIter , TolFeas , TolOpt , Display
%
%   Outputs
%   -------
%   p            : new iterate
%   H            : updated BFGS matrix
%   info struct  : .iter .feas .opt .flag
%
% -------------------------------------------------------------------------
    arguments
        nlp struct
        p   double %column
        H   double
        
        opts.MaxIter (1,1) double = 80
        opts.TolFeas (1,1) double = 1e-8
        opts.TolOpt  (1,1) double = 1e-6
        opts.StepTol  (1,1) double = 1e-2
        opts.Display (1,1) logical = true % false
        opts.CheckGrad (1,1) logical = false
    end

    lb = nlp.lb;  ub = nlp.ub;
    % [f,g]        = nlp.cost(p);          % objective & grad
    % [c , ceq , ~ , dCeq]    = nlp.nonl(p);          % equality rows & Jacobian
    % Aeq = dCeq.';
    % % if it ==1
    %     feas0 = norm(ceq.', inf);
    %     feas = feas0;
    % % end
    for it = 1:opts.MaxIter

        % ----- model evaluation ----------------------------------------
        [f,g]        = nlp.cost(p);          % objective & grad
        [~ , ceq , ~ , dCeq]    = nlp.nonl(p);          % equality rows & Jacobian
        Aeq = dCeq.';
        if it ==1
            feas0 = norm(ceq.', inf);
            feas = feas0;
        end
        if opts.CheckGrad && it == 1            % run once at the first iterate
            disp('Check Gradient:');

            epsFD = 1e-6;
            g_fd  = finiteDiff(@nlp.cost,  p, epsFD);
            A_fd  = finiteDiffJac(@nlp.nonl, p, epsFD);
            
            fprintf('‖g_fd - g‖/‖g_fd‖ = %.2e\n', ...
                    norm(g_fd - g)/max(1,norm(g_fd)));
            fprintf('‖A_fd - Aeq‖/‖A_fd‖ = %.2e\n', ...
                    norm(A_fd - Aeq,'fro')/max(1,norm(A_fd,'fro')));
        end

        % feas = norm(ceq,Inf);
        % optn = norm(g + (Aeq.'*Aeq)\ceq.', Inf);   % KKT residual
        % 
        % if opts.Display
        %     fprintf('#%3d   f=% .3e   feas=% .1e   opt=% .1e\n', ...
        %              it, f, feas, optn);
        % end
        

        % ----- condensed QP search direction ---------------------------
        % d = quadprod(H, g, [],[], Aeq, -ceq, lb-p, ub-p, [])
        %                    % optimoptions('quadprog'));
        % %                    % optimoptions('quadprog','Display','off'));
        [dP,lambda] = solveQP(H,g,Aeq,ceq);
        disp('dp'); disp(norm(dP));
        % feas = norm(ceq.',Inf);
        % % feas = norm(ceq.');
        % optn = norm(g + Aeq.'*lambda , Inf);   % KKT residual
        % % optn = norm(g + Aeq.'*lambda);   % KKT residual
        disp('lambda1');disp(norm(lambda));
        % 
        % if feas < opts.TolFeas && optn < opts.TolOpt
        %     info.flag = 1;  break
        % end
        % if feas > 10*feas0
        %     disp('blow up H');
        %     H = speye(size(H));               % blow-up safeguard
        % end

        % ------------- merit function --------------------------------
        rho  = 10;       % 10                              % penalty weight
        phi0 = f + rho*norm(ceq.',1);
        Ad   = Aeq*dP;                                 % m × 1
        meritSlope = g.'*dP - rho*norm(Ad,1);           % m'(0)
        
        % ------------- backtracking with Armijo ----------------------
        alpha = 1.0;
        % while alpha > 1e-6
        while alpha > 1e-10
            pTry = min( max( p + alpha*dP , lb ) , ub );
            [fTry,gTry]   = nlp.cost(pTry);
            [~ , ceqTry , ~ , AeqTry]     = nlp.nonl(pTry);
            AeqTry = AeqTry.';

            % lambdaTry = -(AeqTry*AeqTry.') \ (AeqTry*gTry);        % m×1 SPD

            phiTry     = fTry + rho*norm(ceqTry.',1);
            dirDer  = g.'*dP - rho*norm(Aeq*dP,1);                 % slope at p
            % dirDer  = gTry.'*dP - rho*norm(AeqTry*dP,1);                 % slope at p

            if phiTry <= phi0 + 1e-4*alpha*dirDer
            % if phiTry <= phi0 + alpha*meritSlope
                break                                    % sufficient decrease
            end
            alpha = alpha*0.5;
            % rho   = rho*2;                               % emphasise feasibility
        end
        
        if alpha <= 1e-10
            info.flag = -2; 
            disp('failed here');
            return 
            % alpha = alpha_prev/2;
        end
        % ---------- BFGS (Powell damped) -----------------------------
        s = pTry - p;
        pOld = p;
        p = pTry;
        % [~,gNew] = nlp.cost(p);
        gNew = gTry;
        Aeq = AeqTry;
        ceq = ceqTry;
        % lambda = lambdaTry;
        % disp('lambda2');disp(norm(lambda));
        % (1) ***Fresh λ consistent with new point*** --------------------
        % lambda = -(Aeq*Aeq.') \ (Aeq*g);            % m×1
        % lambda = -(Aeq*Aeq.') \ (Aeq*gNew);            % m×1
        lambda    = lambda + Aeq*(p-pOld);    % simple least-squares *hint*
        % disp('sdsd'); disp(norm(lambda));
        % (2) ***Correct merit slope uses current Jacobian*** -----------
        % meritSlope = g.'*dP - rho*norm(Aeq*dP,1);   % NOT Aeq_old!
       
        if feas < opts.TolFeas && optn < opts.TolOpt
            info.flag = 1;  break
        end
        y = gNew - g + Aeq.'*lambda;          % full secant
        % y = gNew - g + Aeq.'*lambdaLS;          % full secant

        theta = 0.8 / (1 - s.'*y/(s.'*H*s));
        y     = theta*y + (1-theta)*H*s;      % damping
       
        if s.'*y > 1e-12
            disp('normH'); disp(normest(H));
            Hy = H*s;
            H  = H - (Hy*Hy.')/(s.'*Hy) + (y*y.')/(s.'*y);
        else
            disp('blow up H');

            H  = speye(size(H));              % restart if curvature lost
        end
        f = fTry;
        g = gNew;
        % alpha_prev = alpha;
        optn = norm( gNew + Aeq.'*lambda , Inf );   % ⇽ correct stationarity residual
        % feas = norm( ceq.' , Inf ); 
        feas = norm( ceq.' , 1 ); 
        % lambdaLS = -(Aeq*Aeq.')\(Aeq*gNew);
        % optn     = norm(gNew + Aeq.'*lambdaLS, Inf);
        
        if opts.Display
            fprintf('#%2d  f=%7.3e  feas=%4.1e  opt=%4.1e  α=%4.1e\n', ...
                    it,f,feas,optn,alpha);
        end
        % % ----- projected line-search -----------------------------------
        % alpha = 1.0;
        % 
        % while alpha > 1e-4
        %     pTry = min(max(p + alpha*dP , lb) , ub);  % clip to bounds
        %     [~ , al_ceq , ~ , ~] = nlp.nonl(pTry);
        %     nmSub = norm(al_ceq.',Inf);
        %     % nmSub = norm(al_ceq.');
        %     if nmSub <= (1-opts.StepTol*alpha)*feas
        %         % disp('soved')
        %         break
        %     end
        %     alpha = 0.5*alpha;
        % end
        % if opts.Display
        %     fprintf('#%2d  f=%7.3e  feas=%4.1e  opt=%4.1e  α=%4.1e\n', ...
        %             it,f,feas,optn,alpha);
        % end
        % if alpha <= 1e-4, info.flag = -2; break, end
        % 
        % % ----- BFGS update (full or L-BFGS) ----------------------------
        % s = pTry - p;
        % p = pTry;                                 % accept step
        % [~,gNew] = nlp.cost(p);                   % new gradient
        % y = gNew - g + Aeq.'*lambda;              % secant
        % if s.'*y > 1e-12
        %     Hy = H*s;
        %     H  = H - (Hy*Hy.')/(s.'*Hy) + (y*y.')/(s.'*y);
        % end
    end

    info.iter = it;  info.feas = feas;  info.opt = optn;
    if ~isfield(info,'flag'), info.flag = 0; end
end
function [dP,lambda] = solveQP(H,g,Aeq,ceq)
    % LDLᵀ factorisation -------------------------------------------------
    [L,D,perm] = ldl(H,'vector');
    HinvB  = @(B) HinvTimes(L,D,perm,B);   % works for any RHS
    % HinvB   = @(B) perm(L'\(D\(L\(perm*B))));

    % 1.  H⁻¹ times Aeqᵀ   (n × m)
    HinvAeqT = HinvB( Aeq.' );                % 3624 × 3456
    % HinvAeqT   = Aeq * HinvB( eye(size(H)) );   % (m×n)(n×n) = m×n
    % S = HinvAeqT * Aeq.' ;                           % 3456 × 3456  SPD

    % 2.  Schur complement  S = Aeq * H⁻¹ * Aeqᵀ   (m × m)
    S = Aeq * HinvAeqT;                           % 3456 × 3456  SPD

    % 3.  Condensed RHS     r =  ceq  −  Aeq * H⁻¹ g   (m × 1)
    % r   = ceq.' + Aeq * HinvB(g);               % 3456 × 1
    % r   = ceq.' - Aeq * HinvB(g);               % 3456 × 1
    r   = -(ceq.') - Aeq * HinvB(g);               % 3456 × 1

    % 4.  Solve  S λ = r   (SPD → Cholesky via backslash)
    % lambda = S \ r;
    lambda = -S \ r;
    % lambdaLS = -(Aeq*Aeq.')\(Aeq*g);
                % lambdaTry = -(AeqTry*AeqTry.') \ (AeqTry*gTry);        % m×1 SPD


    % 5.  Recover  Δp = −H⁻¹ ( g + Aᵀ λ )          (n × 1)
    dP     = -HinvB( g + Aeq.'*lambda );
    % dP     = -HinvB( g + Aeq.'*lambdaLS );
end

% function [dP,lambda] = solveQP(H,g,Aeq,ceq)
%     % LDL factorisation with permutation vector
%     % [L,D,perm] = ldl(H,'vector');
%     % 
%     % % helper that works for multiple RHS
%     % HinvB   = @(B) HinvTimes(L,D,perm,B);
%     % 
%     % % AHinv   = Aeq * HinvB( eye(size(H)) );   % (m×n)(n×n) = m×n
%     % % 
%     % % AHinv      = Aeq * HinvB( Aeq.' );           % m×n · n×m
%     % % rhs        = -( ceq.' + AHinv * g );           % ← FIXED SIGN
%     % % lambda     = AHinv \ rhs;                    % SPD solve (Cholesky inside “\”)
%     % % 
%     % % dP         = HinvB( -g - Aeq.'*lambda );     % search direction
%     % % build Schur complement  S = A H⁻¹ Aᵀ
%     % AHinv   = Aeq * HinvB( eye(size(H)) );   % (m×n)(n×n) = m×n
%     % S       = AHinv * Aeq.';                 % m×m  SPD
%     % rhs     = ceq.' - AHinv * g;               % condensed RHS
%     % 
%     % % solve for multipliers and step
%     %  lambda  = S \ rhs;                      % Cholesky inside “\”
%     % dP      = -HinvB( g + Aeq.'*lambda );
% 
%     % % helper that works for multiple RHS
%     % HinvB   = @(B) HinvTimes(L,D,perm,B);
%     % 
%     % % build Schur complement  S = A H⁻¹ Aᵀ
%     % AHinv   = Aeq.' * HinvB( eye(size(H)) );   % (m×n)(n×n) = m×n
%     % S       = AHinv * Aeq;                 % m×m  SPD
%     % rhs     = ceq.' + AHinv * g;               % condensed RHS
%     % 
%     % % solve for multipliers and step
%     % lambda  = -S \ rhs;                      % Cholesky inside “\”
%     % dP      = HinvB( -g - Aeq*lambda );
% end
function X = HinvTimes(L,D,perm,B)
    % solves  H X = B   where  H = P' * L * D * L' * P
    % L, D, perm are from  [L,D,perm] = ldl(H,'vector')
    Y          = B(perm , :);      % apply permutation  P * B
    Y          = L  \ Y;           % forward solve
    Y          = D  \ Y;           % diagonal solve
    Y          = L.' \ Y;           % backward solve
    % Y          = L' \ Y;           % backward solve
    X          = zeros(size(B),'like',B);
    X(perm ,:) = Y;                % X = P' * Y
end

%==================================================================
function g = finiteDiff(fun,x,epsFD)
% g = finiteDiff(@(x)f, x, epsFD)
%   central-difference gradient of a scalar function f(x)

    n = numel(x);   g = zeros(n,1);
    for k = 1:n
        xk  = x;  xk(k) = xk(k) + epsFD;
        f_p = fun(xk);

        xk  = x;  xk(k) = xk(k) - epsFD;
        f_m = fun(xk);

        g(k) = (f_p - f_m) / (2*epsFD);
    end
end
%==================================================================
function A = finiteDiffJac(fun,x,epsFD)
% A = finiteDiffJac(@(x)[c,Jac], x, epsFD)
%   central-difference Jacobian of a vector function c(x)
%   'fun' must return the residual only (no Jacobian).

    [~,c0,~,~] = fun(x);                % vector residual at nominal point
    m  = numel(c0);  n = numel(x);
    A  = zeros(m,n);

    for k = 1:n
        xk       = x;  xk(k) = xk(k) + epsFD;
          
        [~ , ck_plus , ~ , ~]= fun(xk);

        xk       = x;  xk(k) = xk(k) - epsFD;
        [~ , ck_minus , ~ , ~] = fun(xk);

        A(:,k) = (ck_plus - ck_minus) / (2*epsFD);
    end
end
%==================================================================

