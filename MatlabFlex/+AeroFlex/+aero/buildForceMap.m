function FM = buildForceMap(ROM_cts, Na, Nm, disc_method, num_ctrl_surfaces, ctrl_var, gust_opt, force_scale,dataPath)
%BUILDFORCEMAP  Extract D0, D1, B0, B1 … from a Krylov ROM.
%
% Inputs
%   ROM_cts : continuous‑time ss object returned by AeroDiscToCont
%   Na      : number of aero Γ–χ states (cfg.Na)
%
% Outputs
%   FM      : struct with fields
%               D0, D1 , CGamma, Dchi
%               B0, B1 , Bchi
%               AGamma, X

% A = ROM_cts.A;  B = ROM_cts.B;  C = ROM_cts.C;  D = ROM_cts.D;

% ----------------------------------------------------------------------
%   State vector ordering in SHARPy Krylov ROM
%   x = [ Γ ; χ ]    ( Γ length = Na , χ length = Na  → total 2·Na )
%   Inputs  u = [   Γ̇  ;   Γ   ]   (first 2·Na are aero Γ/Γ̇;
%                                    control‑surface inputs follow)
%   Outputs y = modal aerodynamic forces
%
%   So we slice the first 2·Na columns of B and D and partition them:
% ----------------------------------------------------------------------
colsAero = 1 : 2*Na;
colsChi  = 2*Na + (1:3);          % χ̂ rigid‑body (zeros in many ROMs)

A_Gamma = ROM_cts.A;
B_Gamma = ROM_cts.B;
C_Gamma = ROM_cts.C;
D_Gamma = ROM_cts.D;

B0 = B_Gamma(:,1:Nm);
B1 = B_Gamma(:,1+Nm:2*Nm);
D0 = D_Gamma(:,1:Nm);
D1 = D_Gamma(:,1+Nm:2*Nm);

if gust_opt
    % disp('gust enabled');
    Bw = ROM_cts.B(:,2*Nm+1);
    B_delta = ROM_cts.B(:, 2*Nm+2: 2*Nm+1+num_ctrl_surfaces);
    B_ddelta =  ROM_cts.B(:, 2*Nm+2+num_ctrl_surfaces: 2*Nm+1+num_ctrl_surfaces*ctrl_var);

    Dw = ROM_cts.D(:,2*Nm+1);
    D_delta = ROM_cts.D(:, 2*Nm+2: 2*Nm+1+num_ctrl_surfaces);
    D_ddelta =  ROM_cts.D(:, 2*Nm+2+num_ctrl_surfaces: 2*Nm+1+num_ctrl_surfaces*ctrl_var);


    % check indexing:
    B_err = norm(ROM_cts.B - [B0 B1 Bw B_delta B_ddelta]);
    D_err = norm(ROM_cts.D - [D0 D1 Dw D_delta D_ddelta]);

else
    % disp('gust disabled');

    Bw = zeros(Na,1);
    Dw = zeros(Nm,1);
    B_delta = ROM_cts.B(:, 2*Nm+1: 2*Nm+num_ctrl_surfaces);
    B_ddelta =  ROM_cts.B(:, 2*Nm+1+num_ctrl_surfaces: 2*Nm+num_ctrl_surfaces*ctrl_var);
    
    D_delta = ROM_cts.D(:, 2*Nm+1: 2*Nm+num_ctrl_surfaces);
    D_ddelta =  ROM_cts.D(:, 2*Nm+1+num_ctrl_surfaces: 2*Nm+num_ctrl_surfaces*ctrl_var);
    
    % check indexing:
    B_err = 0; %norm(SS_sub.B - [B0 B1 B_delta B_ddelta]);
    D_err = 0; %norm(SS_sub.D - [D0 D1 D_delta D_ddelta]);
end   

if B_err ~= 0 || D_err ~=0
    error('Indexing Not Correct For B or D')
end
% D0 (stiffness)     corresponds to Γ component of input
% D1 (damping)       corresponds to Γ̇ component
FM.D0 = D0/force_scale;
FM.D1 = D1/force_scale;

FM.B0 = B0;
FM.B1 = B1; 

% % χ̂ placeholders (zeros if the ROM has no rigid‑body states)
% FM.Dchi = zeros(size(FM.D0,1),3);
% FM.Bchi = zeros(size(B_Gamma,1),3);

% CGamma and AGamma are the aero output / state blocks
FM.C_Gamma = C_Gamma/force_scale;
FM.A_Gamma = A_Gamma;

% Gust and Control Partition
FM.Bw = Bw;
FM.B_delta = B_delta;
FM.B_ddelta = B_ddelta;

FM.Dw = Dw/force_scale;
FM.D_delta = D_delta/force_scale;
FM.D_ddelta = D_ddelta/force_scale;

% Gravity / rigid‑body coupling (zero for now) Temp 
% FM.X = zeros(3,size(C_Gamma,1));


% OrientFile = fullfile(dataPath, 'OrientMats.mat');
% OrientStruct =  load(OrientFile); X = OrientStruct.X;
% FM.X = X;



end
