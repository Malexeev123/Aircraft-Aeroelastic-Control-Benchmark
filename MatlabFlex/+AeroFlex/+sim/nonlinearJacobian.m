function [Nq,Nu] = nonlinearJacobian(x_n, idx, ...
                                     par)
%NONLINEARJACOBIAN  Analytic Jacobians ∂N/∂q  &  ∂N/∂u  for Artola ROM
%
%   [Nq,Nu] = nonlinearJacobian(q1,q2,qxi,u, ...
%                               G1,G2,Gg,Gxi, ...
%                               Bdelta,Bgust, Naero)
%
%   q1,q2    : Nm×1          
%   qxi      : (Nm+1)×1      rigid‑body ξ coordinates
%   u        : m×1           control/disturbance  [ δ̇ ; w_g ]
%
%   G1,G2    : Nm×Nm×Nm      Γ₁, Γ₂ tensors
%   Gg       : Nm×(Nm+1)×(Nm+1)   Γ_g tensor
%   Gxi      : (Nm+1)×Nm×(Nm+1)   Γ_ξ tensor
%   Bdelta   : n×1           state‑forcing vector for flap rate
%   Bgust    : n×1           state‑forcing vector for gust
%   Naero    : integer       aerodynamic states (≥0)
%
%   Nq       : n×n           Jacobian w.r.t. state
%   Nu       : n×m           Jacobian w.r.t. input   (here m = 2)
%  ────────────────────────────────────────────────────────────────────
q1      = x_n(idx.q1 );            % Nm
q2      = x_n(idx.q2 );            % Nm
qxi    = x_n(idx.qxi);
G1 = par.Gamma1;
G2 = par.Gamma2;
Gg = par.Gamma_g;
Gxi = par.Gamma_xi;

Naero = par.Na;

Nm   = numel(q1);
Nxi  = Nm+1;
nx   = 3*Nm + 1 + Naero+ 3;              % total state dimension
% m    = 2;                             % # inputs  (δ̇ , w_g)
ctrlIdx1 = size(par.Bdel, 2);
ctrlIdx2 = size(par.Bddel,2);
m    = size(par.Bw,2) + ctrlIdx1 + ctrlIdx2;                             % # inputs  (δ̇ , w_g)
% Nq = zeros(x_n,x_n);
% ---------------------------------------------------------------------
% Helper to contract Γ‑tensors  Σ_l q(l) Γ^{l}_{⋅⋅}
% ---------------------------------------------------------------------
% contract = @(T, q) squeeze( pagemtimes( ...
%            permute(T,[1 3 2]), ...    % (i,l,k) or (i,l,j)
%            repmat(reshape(q,[],1,1),1,1,size(T,3)) ) );
% contract = @(T,q) ...
%     squeeze( pagemtimes( permute(T,[1 3 2]), ...
%                          repmat(reshape(q,[],1,1), 1,1,size(permute(T,[1 3 2]),3)) ) );
% contract = @(T,q) ...
%     squeeze( pagemtimes(permute(repmat(reshape(q,1,[],1), 1,1,size(permute(T,[2 3 1]),3)),[1 3 2] ), ...
%                          T ) );
% contract2 = @(T,q) ...
%     squeeze( pagemtimes(permute(repmat(reshape(q,[],1,1), 1,1,size(permute(T,[2 3 1]),3)),[2 3 1] ), ...
%                          T ) );

% ---------- ∂N/∂q1  (Γ₁ block) ---------------------------------------
% % termA = contract(G1 , q1);            % Σ_l q1(l) Γ^{l}_{ik}
% termA = permute(pagemtimes(permute(repmat(q1,1,1,size(G1,1 )),[2 1 3]), permute(G1,[1 3 2])),[2 3 1]);            % Σ_l q1(l) Γ^{l}_{ik}
% % termA = contract(G1 , reshape(q1,1,1,[]));            % Σ_l q1(l) Γ^{l}_{ik}
% % termB = pagemtimes(G1 , reshape(q1,1,1,[]));   % Σ_j q1(j) Γ^{k}_{ij}
% termB = permute(pagemtimes(permute(repmat(q1,1,1,size(G1,1 )),[2 1 3]), permute(G1,[2 1 3])),[2 3 1]);
% % termB = pagemtimes(G1 , reshape(q1,[],1,1));   % Σ_j q1(j) Γ^{k}_{ij}
% termB = squeeze(termB);
% Nq11  = -( termA + termB );           % Nm × Nm
% 
% % ---------- ∂N/∂q2  (Γ₂ block) ---------------------------------------
% termC = contract(G2 , q2);
% % termD = pagemtimes(G2 , reshape(q2,1,1,[]));
% termD = pagemtimes(G2 , reshape(q2,[],1,1));
% termD = squeeze(termD);
% Nq12  = -( termC + termD );           % Nm × Nm
% 
% % ---------- ∂N/∂qξ (Γ_g block) --------------------------------------
% termE = contract(permute(Gg,[3 1 2]) , qxi);           % Nm×(Nm+1)
% termF = contract(permute(Gg,[2 1 3]), qxi);
% Nq13  =   ( termE + termF );          % Nm × (Nm+1)
% 
% % ---------- Γ₂ᵀ block  ∂N/∂q2  in 2nd rows ---------------------------
% % termG = contract(permute(G2,[3 2 1]), q2);
% termG = contract(permute(G2,[2 3 1]), q2);
% % termH = pagemtimes(permute(G2,[3 2 1]), reshape(q1,1,1,[]));
% termH = pagemtimes(permute(G2,[3 2 1]), reshape(q1,[],1,1));
% termH = squeeze(termH);
% Nq21  = termG ;           % Nm × Nm
% Nq22 = termH ;                    % Nm × Nm
% 
% % ---------- Γ_ξ block  ∂N/∂qξ  in 3rd rows ---------------------------
% % termI = contract(Gxi , reshape(q1,1,[],1));           % (Nm+1) × (Nm+1)
% % termI = contract(permute(Gxi,[1 2 3]) , reshape(q1,[],1,1));           % (Nm+1) × (Nm+1)
% % termI = contract(Gxi , reshape(q1,[],1,1));           % (Nm+1) × (Nm+1)
% % termI = contract(permute(Gxi,[1 3 2]) , reshape(q1,1,[],1));           % (Nm+1) × (Nm+1)
% termI = permute(pagemtimes(permute(repmat(q1,1,1,size(Gxi,1 )),[2 1 3]), permute(Gxi,[2 1 3])),[2 3 1]);           % (Nm+1) × (Nm+1)
% % termI = permute(pagemtimes(permute(repmat(q1,1,1,size(Gxi,1 )),[2 1 3]), permute(Gxi,[2 1 3])),[3 2 1]);           % (Nm+1) × (Nm+1)
% 
% % termI = contract(Gxi , q1);           % (Nm+1) × (Nm+1)
% % termJ = contract(permute(Gxi,[1 3 2]), qxi);
% % termJ = contract(Gxi, qxi);
% termJ = permute(pagemtimes((permute(repmat(qxi,1,1,size(Gxi,2 )),[2 1 3])), permute(Gxi, [1 3 2])),[2 3 1]);
% Nq31  =   termJ;          % (Nm+1)×(Nm+1)
% Nq33 = termI;
%% Attempt Manual

% termA2 = zeros(size(G1,1), size(G1,2));
% termB2 = [];
% for i = 1: size(G1,3)
%     sumA2 = q1(i)*G1(:,:,i);
%     termA2 = termA2 +sumA2;
%     mult = G1(:,:,i)* q1;
%     termB2 = [termB2, mult];
% end
%%
% % Γ₁  (20×20×20) ---------------------------------------------
termA = gammaContract( permute(G1,[1 3 2]), q1 );  % Σ_l q1(l) Γ^l_{ik} Actual
termB = gammaContract( G1 , q1 );                  % Σ_j q1(j) Γ^k_{ij}

% Nq11  = -(termA + termB);
Nq11  = -termA - termB;
% Nq11  = -termA2 - termB2;

% Γ₂  (20×20×20)  same pattern -------------------------------
%% Attempt Manual
% termC2 = squeeze(pagemtimes(G2, 'none', q2, 'none'));
% termC2 = zeros(size(G2,1), size(G2,2));
% termD2 = [];
% for i = 1: size(G2,3)
%     sumC2 = q2(i)*G2(:,:,i);
%     termC2 = termC2 +sumC2;
%     mult = G2(:,:,i)* q2;
%     termD2 = [termD2, mult];
% end

termC = gammaContract( permute(G2,[1 3 2]), q2 ); % Actual
termD = gammaContract( G2 , q2 );

% Nq12  = -(termC + termD);
Nq12  = - termC - termD;
% Nq12  = -(termC2 + termD2);
% Nq12  = -termC2 -termD2;

%% Γ_g (20×21×21)  qs have length 21 --------------------------
% Attempt Manual
% termE2 = zeros(size(Gg,1), size(Gg,2));
% termF2 = [];
% for i = 1: size(Gg,3)
%     sumE2 = qxi(i)*Gg(:,:,i);
%     termE2 = termE2 +sumE2;
%     mult = Gg(:,:,i)* qxi;
%     termF2 = [termF2, mult];
% end

termF = gammaContract( Gg , qxi );
termE = gammaContract( permute(Gg,[1 3 2]), qxi );


Nq13  =  (termE + termF);
% Nq13  =  (termE2 + termF2);

% Γ_ξ (21×20×21)  second dim = 20, permuted second = 21 -------

%% Attempt Manual
% termI2 =  squeeze(pagemtimes(Gxi, 'none', q1, 'none'));
% termI2 =  squeeze(pagemtimes(Gxi, 'none', qxi, 'none'));
% termI2 = zeros(size(Gxi,1), size(Gxi,2));
% termJ2 = [];
% for i = 1: size(Gxi,3)
%     sumI2 = qxi(i)*Gxi(:,:,i);
%     termI2 = termI2 +sumI2;
%     mult = Gxi(:,:,i)* q1;
%     termJ2 = [termJ2, mult];
% end
%% LEGACY 2025
% termJ = gammaContract( Gxi , q1  );          % Σ_j q1(j) Γ^k_{ij} % Actual
% termI = gammaContract( permute(Gxi,[1 3 2]), qxi ); % Σ_l qξ(l) Γ^l_{ik}
% Nq31  =  termI;
% Nq33 = termJ;
%% NEW 2026 JAN
termJ = gammaContract( Gxi , qxi  );          % Σ_j q1(j) Γ^k_{ij} % Actual
termI = gammaContract( permute(Gxi,[1 3 2]), q1 ); % Σ_l qξ(l) Γ^l_{ik}
Nq31  =  termJ;
Nq33 = termI;

% Nq33  =  termJ2;
% Nq31 = termI2;

%% Attempt Manual
% termG2 = zeros(size(G2,1), size(G2,2));
% termH2 = [];
% for i = 1: size(G2,3)
%     sumG2 = q2(i)*(G2(:,:,i).');
%     termG2 = termG2 +sumG2;
%     mult = (G2(:,:,i).')* q1;
%     termH2 = [termH2, mult];
% end

G2T2 = permute(G2,[2 3 1]);                    % (k , j , l) → (i , l , k)
termG = gammaContract( G2T2 , q2 );            % Nm × Nm   (i rows)
% ----- termH :  Σ_j q1(j) Γ₂^{Tj}_{ik}  ---------------------------
termH = gammaContract( pagetranspose(G2) , q1 );
% ----- assemble block --------------------------------------------
Nq21 = termG ;                    % Nm × Nm
Nq22 = termH ;                    % Nm × Nm
% Nq21 = termG2 ;                    % Nm × Nm
% Nq22 = termH2 ; 
% ---------- Assemble structural Jacobian -----------------------------
NqS = zeros(numel(x_n),numel(x_n) );                  % structural part

% rows 1:Nm
NqS(1:Nm , 1:Nm           ) = Nq11;
NqS(1:Nm , Nm+1:2*Nm      ) = Nq12;
NqS(1:Nm , 2*Nm+1:3*Nm+1  ) = Nq13;

% rows Nm+1:2*Nm
NqS(Nm+1:2*Nm , 1:Nm      ) = Nq21;
NqS(Nm+1:2*Nm , Nm+1:2*Nm       ) = Nq22;

% (Γ₂ᵀ off‑diagonal blocks omitted here if zero)

% rows 2*Nm+1:3*Nm+1
NqS(2*Nm+1:3*Nm+1 , 1:Nm ) = Nq31;
NqS(2*Nm+1:3*Nm+1 , 2*Nm+1:3*Nm+1 ) = Nq33;


% ---------- Append aerodynamic zero blocks ---------------------------
% Nq = blkdiag(NqS, zeros(Naero));
Nq = NqS;
% ---------- Input Jacobian (linear maps) ------------------------------
Nu = zeros(nx, nx+m);
% Nu(:,1) = Bdelta;      % flap‑rate
% Nu(:,2) = Bgust;       % gust
Nu(idx.q1, nx+1) = par.gustSet.a*par.gustSet.t_inf * par.Dw;
Nu(idx.qGam, nx+1) = par.Bw;

Nu(idx.q1,nx+ 2:nx+1 + ctrlIdx1) = par.gustSet.a*par.Ddel;
Nu(idx.q1, nx+2 + ctrlIdx1:nx+1 + ctrlIdx1 + ctrlIdx2) = par.scaleAero*par.Dddel;

Nu(idx.qGam, nx+2:nx+1 + ctrlIdx1) = (1/par.t_inf)*par.Bdel;
Nu(idx.qGam, nx+2 + ctrlIdx1:nx+1 + ctrlIdx1 + ctrlIdx2) =  par.Bddel;

end



function M = gammaContract(T,q)
    %GAMMACONTRACT   M = Σ_l q(l) · T(:,l,:)  (matrix output p×k)
    %
    %   T : p × l × k        (second dimension is summation index)
    %   q : l × 1
    %   M : p × k
    %
    lDim = size(T,2);
    if numel(q) ~= lDim
        error('gammaContract: q length %d ≠ 2nd dim of T (%d).',...
               numel(q), lDim);
    end
    % replicate q across pages  → l × 1 × k
    Qpage = repmat( reshape(q,lDim,1,1) , 1,1,size(T,3) );
    % page‑wise multiply: (p×l×k) * (l×1×k) = (p×1×k)
    M3 = pagemtimes(T , Qpage);
    % squeeze singleton 2nd dim, reshape to p×k
    M  = reshape(M3, size(T,1), size(T,3));
end