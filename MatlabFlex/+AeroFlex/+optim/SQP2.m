function SQP2
    % ------------------------------------------------------------
    %  SQP loop   (Gauss–Newton + BFGS)   replaces fmincon
    % ------------------------------------------------------------
    maxIter = 80;                    % <-- tweak as needed
    tolOpt  = 1e-6;                  % gradient tolerance
    tolFeas = 1e-8;                  % equality tolerance
    
    p  = obj.z0;                     % warm start
    H  = speye(numel(p));            % initial Hessian
    
    % feas = inf;
    % optn = inf;
    for it = 1:maxIter
    
        % ---- evaluate cost, grad, constraints ------------------
        [J , gradJ]          = nlp.cost(p);      % analytic grad
        [c , ceq , ~ , dCeq] = nlp.nonl(p);      % ceq, gradceq (Jac)
        

        % if feas < tolFeas && optn < tolOpt,   break,  end

    
        % ---- QP search direction  ------------------------------
        [deltaP,lambda] = solveQP(H, gradJ, dCeq, ceq);
        feas = norm(ceq,Inf);
        % optn = norm(gradJ + dCeq.'*jacobiVec(ceq), Inf); % first-order cond
        optn          = norm(g + dCeq.'*lambda , Inf);

        if feas < tolFeas && optn < tolOpt,   break,  end

        % ---- line search with bound clipping -------------------
        alpha = 1.0;
        while alpha > 1e-3
            pNew = min(max(p + alpha*deltaP, nlp.lb), nlp.ub);
            ceqNew = nlp.nonl(pNew);
            if norm(ceqNew,Inf) < 0.8*feas,  break, end   % Armijo on feas
            alpha = 0.5*alpha;
        end
        if alpha <= 1e-3, warning('SQP-line-search failed');  break, end
    
        % ---- BFGS Hessian update -------------------------------
        s = pNew - p;                         % step
        y = (nlp.costGrad(pNew) - gradJ) ...  % new grad minus old
            + dCeq.'*(lambda - jacobiVec(ceq));
        if s.'*y > 1e-12
            Hy = H*s;
            H  = H - (Hy*Hy.')/(s.'*Hy) + (y*y.')/(s.'*y);
        end
    
        % advance
        p = pNew;
    end
    
    if it == maxIter
        warning('nMHE:SQP','SQP hit maxIter, keep last iterate');
    end

    % p    = obj.z0;                               % warm-start
    % H    = speye(numel(p));                      % initial Hessian  (can recycle)
    % 
    % % ---------- parameter -------------------------------------------------
    % maxIter = 120;                               % typical conv. 40–80
    % tolFeas = 1e-8;   tolOpt = 1e-6;
    % 
    % for it = 1:maxIter
    %     % ---- evaluate model ---------------------------------------------
    %     [J , g ]              = nlp.cost(p);         % gradJ  (analytic)
    %     [ceq , Aeq ]          = nlp.nonl(p);         % ceq , gradceq
    % 
    %     feas = norm(ceq,Inf);
    %     optn = norm(g + Aeq.'*(Aeq\ceq) , Inf);      % F.O. condition
    % 
    %     if feas < tolFeas && optn < tolOpt, break, end
    % 
    %     % ---- condensed QP  Δp = argmin ½ΔpᵀHΔp+gᵀΔp  s.t. AeqΔp=-ceq ----
    %     [dP,lambda] = solveQP(H,g,Aeq,ceq);          % helper shown next
        % 
        % % ---- project to box & line search (‘funnel’) --------------------
        % alpha = 1.0;
        % while alpha > 1e-3
        %     pTry = min(max(p + alpha*dP , nlp.lb), nlp.ub);   % clip
        %     if norm(nlp.nonl(pTry),Inf) <= (1-1e-2*alpha)*feas
        %         break
        %     end
        %     alpha = 0.5*alpha;
        % end
        % if alpha <= 1e-3, warning("SQP stalled"); break, end
        % 
        % % ---- BFGS / L-BFGS update  (curvature test) ---------------------
        % s = pTry - p;          p = pTry;
        % [~, gNew]  = nlp.cost(p);
        % y = gNew - g + Aeq.'*lambda;
        % if s.'*y > 1e-10
        %     Hy = H*s;
        %     H  = H - (Hy*Hy.')/(s.'*Hy) + (y*y.')/(s.'*y);    % full BFGS
            % >>> for L-BFGS store (s,y) pairs instead  <<<
    %     end
    % end
    
    % ---------- extract estimates & warm-start next call ------------------
    obj.z0   = p;                       % store full vector
    obj.xhat = p( obj.idx.x{obj.Ne+1} );% newest node  (=k)
    obj.what = p( [obj.idx.w{:}] );     % 24-element history
    xhat     = obj.xhat;
end

function [dP,lambda] = solveQP(H,g,A,r)
    % condensed Schur complement solve  (Gill & Saunders SC-QP 2005)
    [L,D,perm] = ldl(H,'vector');
    Hinvtimes  = @(b) perm*(L'\(D\(L\(perm.'*b))));   % H⁻¹ b
    AHinv      = A * Hinvtimes(speye(size(H)));
    S          = AHinv * A.';                 % SPD  (Schur)
    rhs        = r + AHinv * g;
    lambda     = -S \ rhs;                    % Cholesky implicitly
    dP         = Hinvtimes(-g - A.'*lambda);
end