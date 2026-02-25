function flutterAnalysis(beam,aero,base,cfg)
%FLUTTERANALYSIS  V‑g and V‑f curves for the current ROM.
%
% cfg fields:
%   Umin,Umax,dU : sweep range [m/s]
%   rho          : air density [kg/m³]
%   savePlots    : true/false
%
% OUTPUTS: plots in current figure and optional PNG file

assert(all(isfield(cfg,{'Umin','Umax','dU','rho'})),'cfg incomplete');

Nm  = beam.Nm;                              % # structural modes
Ns  = size(aero.A_Gamma,1);                 % # aerodynamic states
Nx  = 2*Nm + Ns;                            % total states (no quaternions)

% pre‑allocate storage
nU  = numel(cfg.Umin:cfg.dU:cfg.Umax);
eigMat = complex(zeros(Nx,nU));

% structural & aero matrices (modal coordinates)
[M,K,Sigma] = beam.linearModelMatrices;     % M, K, damping Σ
[Aa,Ba,Ca]  = aero.linearModelMatrices;     % UVLM ROM

% coupling matrices (from base)
[Bq0,Bq1] = base.AICstruct;                 % modal displ/vel to aero

%% sweep loop ------------------------------------------------------------
idx = 0;
for U  = cfg.Umin:cfg.dU:cfg.Umax
    idx = idx+1;

    % dynamic‑pressure‑scaled aero matrices  q = ½ρU²
    qInf = 0.5*cfg.rho*U^2;
    A11  = zeros(Nm);                       % q0′ = q1
    A12  = eye(Nm);
    A21  = -K - Sigma*qInf*eye(Nm);         % simple follower damping
    A22  = zeros(Nm);

    % coupling with aerodynamic states
    A31 = zeros(Ns,2*Nm);
    A31(:,1:Nm) = Ba.q0;                    % qΓ′ = AΓ qΓ + B q0 + B q1
    A31(:,Nm+1:end) = (Ba.q1*qInf);

    A = [ A11  A12          zeros(Nm,Ns) ;
          A21  A22          -Ca.q0*qInf ;
          A31                Aa ];

    eigMat(:,idx) = eig(A);
end

%% post‑process: V‑g and V‑f ---------------------------------------------
V = cfg.Umin:cfg.dU:cfg.Umax;
wn = imag(eigMat);   % rad/s
z  = -real(eigMat)./abs(eigMat);   % damping ratio

figure('Color','w','Position',[100 100 800 600]); hold on
for k = 1:size(eigMat,1)
    plot(V,z(k,:),'b');
end
xlabel('U_\infty  [m/s]'); ylabel('Damping ratio  ζ');
title('V‑g curve'); grid on
yline(0,'k-','Flutter boundary');

figure('Color','w','Position',[950 100 800 600]); hold on
for k = 1:size(eigMat,1)
    plot(V,wn(k,:)/(2*pi),'b');
end
xlabel('U_\infty  [m/s]'); ylabel('Frequency  [Hz]');
title('V‑f curve'); grid on

if isfield(cfg,'savePlots') && cfg.savePlots
    exportgraphics(gcf,'Flutter_Vf.png','Resolution',300);
    figure(1), exportgraphics(gcf,'Flutter_Vg.png','Resolution',300);
    disp('Flutter plots saved.');
end
end
