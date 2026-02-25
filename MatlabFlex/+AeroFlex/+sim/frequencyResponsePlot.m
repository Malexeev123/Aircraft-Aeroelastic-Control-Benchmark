function frequencyResponsePlot(beam,aero,base,cfg)
%FREQUENCYRESPONSEPLOT  Bode‑style FRF for tip or intermediate output.
%
% cfg fields:
%   U            : operating airspeed  [m/s]
%   wMin,wMax    : frequency sweep     [rad/s] or [Hz] (auto)
%   nPts         : # points
%   outNode      : 'tip' or node index  (vertical displacement)
%   inType       : 'gust'  (vertical LE gust)
%                  'flap'  (flap deflection rate)
%
% NOTE: assumes linearised ROM around undeformed config.

% linearised structural matrices
[M,K,Sigma] = beam.linearModelMatrices;
Nm = beam.Nm;

% aerodynamic ROM matrices at this U
qInf  = 0.5*cfg.rho*cfg.U^2;
[Aa,Ba,Ca] = aero.linearModelMatrices;
Ns = size(Aa,1);

% assemble state‑space   x = [q0;q1;qΓ]
A = [ zeros(Nm)   eye(Nm)            zeros(Nm,Ns);
     -K          zeros(Nm)          -Ca.q0*qInf ;
      Ba.q0      Ba.q1*qInf          Aa          ];

if strcmpi(cfg.inType,'gust')
    B = [ zeros(2*Nm,1) ; Ba.w * qInf ];
    inLabel = 'vertical gust  w_g [m/s]';
elseif strcmpi(cfg.inType,'flap')
    B = [ zeros(2*Nm,1) ; Ba.dDel * qInf ];
    inLabel = 'flap rate  \dot\delta [rad/s]';
else
    error('Unknown inType');
end

% output matrix: tip vertical displacement y = e_tip^T * Phi0 * q0
if ischar(cfg.outNode) && strcmpi(cfg.outNode,'tip')
    node = size(beam.fem.coordinates,1);     % last node
else
    node = cfg.outNode;
end
rows = AeroFlex.core.localDofIndex(beam.fem.keepDofs,(node-1)*6+2); % global Y disp row
C = zeros(1,2*Nm+Ns);
C(1,1:Nm) = beam.phi0(rows,:);              % selects q0
D = zeros(1,1);

sys = ss(A,B,C,D);

% frequency vector
if cfg.wMax>50       % assume rad/s if >50
    w = logspace(log10(cfg.wMin),log10(cfg.wMax),cfg.nPts);
    [mag,phs] = bode(sys,w);
    f = w/(2*pi);    % convert to Hz for plotting
else
    f = logspace(log10(cfg.wMin),log10(cfg.wMax),cfg.nPts);
    [mag,phs] = bode(sys,2*pi*f);
    mag = squeeze(mag); phs = squeeze(phs);
end

figure('Color','w','Position',[200 200 900 700])
subplot(2,1,1)
semilogx(f,20*log10(mag),'b','LineWidth',1.2)
ylabel('|H(j\omega)|  [dB]'); grid on
title(sprintf('Frequency response  (U = %.1f m/s, output = %s)',...
      cfg.U,cfg.outNode))
subplot(2,1,2)
semilogx(f,phs,'r','LineWidth',1.2)
ylabel('Phase  [deg]'); xlabel('Frequency  [Hz]'); grid on

if isfield(cfg,'savePlots') && cfg.savePlots
    savefig(sprintf('FRF_%s_U%.0f.fig',cfg.outNode,cfg.U));
end
end
