function result = smoketest_uniform_cantilever()
% Build a uniform Euler–Bernoulli cantilever, add small Kelvin–Voigt,
% excite first bending, measure log decrement and frequency.

% Geometry & properties
L   = 1.0;     % m
Nn  = 21;      % nodes
A   = 2.0e-4;  % m^2
Iy  = 3.0e-8;  % m^4
Iz  = 3.0e-8;  % m^4
GJ  = 1.0;     % (unused EB) keep numeric
E   = 70e9;    % Pa
rho = 2700;    % kg/m^3
eta_f = 5e3;   % N*s/m^2  (tiny)
eta_m = 5e-3;  % N*m*s    (tiny)

% 1) build cantilever beam mesh
[fem, keepDofs, chain, elemLen, Eall] = AeroFlex.beam.make_uniform_beam(L,Nn,'cantilever',E,A,Iy,Iz,GJ,rho);
fem.E = Eall;
% 2) linear M,K, modes
[M,K,Vr,w] = AeroFlex.beam.linear_modal_analysis(fem, keepDofs); % implement thin wrapper to your existing build
Mred = Vr.'*M*Vr; scl = 1./sqrt(diag(Mred)); phi0_mn = Vr.*scl.'; omega = w(:);

% 3) damping
opts = struct('eta_f',eta_f,'eta_m',eta_m);
Dlin = AeroFlex.beam.build_linearized_kv_damping(fem, keepDofs, chain, elemLen, Eall, opts);
Dr   = AeroFlex.beam.project_kv_damping_to_modes(phi0_mn, Dlin);

% 4) simulate single-mode free response: q̈ + 2ζω q̇ + ω² q ≈ 0
j=1;
zeta_meas = Dr(j,j)/(2*omega(j)); % since Mhat=I
dt=1e-4; T=0.5; t=(0:dt:T).';
q = zeros(length(t),1); qd=q; q(1)=1e-4;  % small IC
for k=2:length(t)
    % simple Newmark/central difference on scalar SDOF with viscous c=2ζω
    c = 2*zeta_meas*omega(j);
    a = ( -c*qd(k-1) - omega(j)^2*q(k-1) );
    qd(k) = qd(k-1) + dt*a;
    q(k)  = q(k-1) + dt*qd(k);
end

% estimate ζ from peaks
[pks,locs] = findpeaks(abs(q), 'MinPeakDistance', round(0.5/(omega(j)/2/pi)/dt));
if numel(pks)>=2
    delta = log(pks(1)/pks(2));
    zeta_from_logdec = delta/sqrt(4*pi^2 + delta^2);
else
    zeta_from_logdec = NaN;
end

result = struct('omega1',omega(1), 'zeta_from_Dr', zeta_meas, 'zeta_from_logdec', zeta_from_logdec);
disp(result)
end
