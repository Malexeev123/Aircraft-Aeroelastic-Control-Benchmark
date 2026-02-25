
function params = paramsRigid_PazyUAV()
% paramsRigid_PazyUAV 
% Defines the rigid-body masses, inertias, and locations for the Pazy wing UAV configuration.
% Returns a structure with fields for fuselage, tail, battery, and tip masses, including their 
% mass values, inertia tensors, and position vectors relative to the wing root (elastic axis origin).
%
% All units are SI (meters, kilograms, kg-m^2).
    
    % Fuselage parameters
    params.fuselage.mass = 0.20;  % mass of fuselage (kg) – estimate
    fus_length = 0.40;           % fuselage length (m)
    fus_radius = 0.025;          % fuselage radius (m) (cylinder radius ~0.025 m => 5 cm diameter)
    % Fuselage center of mass: assume fuselage is uniformly distributed and wing attaches at fuselage center
    params.fuselage.pos = [0, 0, 0];  % (x,y,z) of fuselage CG relative to wing root (here origin at wing-fuselage attach)
    % Fuselage inertia (approximated as a uniform cylinder about its center, aligned along body X-axis)
    M = params.fuselage.mass;
    R = fus_radius;
    L = fus_length;
    Ixx = 0.5 * M * R^2;                  % moment about fuselage long axis (roll)
    Iyy = (1/12) * M * (3*R^2 + L^2);     % moment about transverse axis (pitch)
    Izz = Iyy;                            % moment about vertical axis (yaw), same for symmetric cylinder
    params.fuselage.inertia = diag([Ixx, Iyy, Izz]);  % inertia matrix (3x3) in body axes (aligned with global axes here)
    
    % Tail mass (horizontal tail + any tail boom mass not included in fuselage)
    params.tail.mass = 0.03;  % kg (estimate for tail surfaces)
    params.tail.j = .01;
    % Assume tail CG is located at fuselage end (tail) with no vertical offset
    params.tail.pos = [-0.20, 0, 0];  % at tail, 0.20 m behind wing root (negative X since wing root is mid fuselage and positive X is forward)
    % Tail inertia – treat as point mass (negligible rotational inertia about its CG for this model)
    % params.tail.inertia = diag([0, 0, 0]);
    params.tail.inertia = diag([params.tail.j, .5*params.tail.j, .5*params.tail.j]);
    
    % Battery mass
    params.battery.mass = 0.10;  % kg (approximate battery weight)
    % Assume battery is located towards the nose inside fuselage
    params.battery.pos = [0.10, 0, 0];  % 0.10 m forward of wing root along fuselage (inside the nose section)
    params.battery.inertia = diag([0, 0, 0]);  % treat as point mass for simplicity
    
    % Wing tip masses (one on each wing tip)
    params.tip.mass = 0.025;  % kg for each tip mass (small weights attached at wing tips)
    % Positions of tip masses relative to wing root.
    % Wing span is 0.55 m (half-span 0.275 m). Tip masses are located at the wing quarter-chord at each tip.
    half_span = 0.275;
    % Assuming no sweep or dihedral, wing tip quarter-chord lies roughly at same X,Z as root (x=0, z=0 in undeformed state).
    params.tip.pos = [0, -half_span, 0;    % left tip (negative Y direction)
                      0,  half_span, 0];   % right tip (positive Y direction)
    % Tip mass inertia (point masses at tip, negligible rotational inertia)
    params.tip.inertia = diag([0, 0, 0]);
    
    % (Optional) Fuselage geometry for visualization purposes
    params.fuselageGeom.length = fus_length;
    params.fuselageGeom.radius = fus_radius;
    
    % ----- user-tunable geometric helpers ---------------------------------
    c      = 0.10;             % wing chord [m]
    ea     = 0.42*c;           % elastic-axis x-offset from LE
    bHalf  = 0.55;             % half-span to tip mass
    xWingRoot = 0.0;           % place root at X = 0
    zWing     = 0.0;           % mid-plane
    
    % ======================================================================
    %  1. Lumps
    % ======================================================================
    
    % Fuselage  -------------------------------------------------------------
    fus.m  = 0.60;
    fus.r  = [ xWingRoot-ea/2 ,   0 , zWing ];      % put CM slightly ahead of EA
    fus.J  = diag([7.1e-3  , 2.3e-2 , 2.3e-2]);
    
    % Battery (inside fuse) -------------------------------------------------
    bat.m  = 0.18;
    bat.r  = fus.r + [0.08 , 0 , 0];                % 8 cm forward
    bat.J  = diag([3.3e-4 , 1.9e-3 , 2.2e-3]);
    
    % % Horizontal tail  ------------------------------------------------------
    % tail.m = 0.12;
    % tail.r = fus.r + [-0.60 , 0 , 0];               % 0.6 m aft
    % tail.J = diag([2.0e-3 , 1.1e-3 , 3.0e-3]);
    % 

    % values if available.
    tail.m = 0.06;                                      % mass (kg)
    tail.r = fus.r + [-0.60, 0, 0];                     % 0.6 m aft (m)
    tail.J = diag([1.0e-3, 5.0e-4, 1.5e-3]);             % inertia (kg·m²)

    % Vertical tail (fin): a small vertical stabiliser contributes to
    % directional stability.  We model it as a lumped mass located aft
    % and slightly above the wing plane.  The inertia tensor is
    % approximate and assumes the fin is slender.
    fin.m  = 0.03;                                      % mass (kg)
    fin.r  = fus.r + [-0.60, 0, 0.05];                  % 0.6 m aft, 5 cm up (m)
    fin.J  = diag([5.0e-4, 1.0e-3, 1.0e-3]);             % inertia (kg·m²)

    % Wing tip masses (actuators) ------------------------------------------
    tipR.m = 0.035;  tipL.m = tipR.m;
    tipR.r = [ xWingRoot + bHalf*tan(0) ,  +bHalf , zWing ];
    tipL.r = [ xWingRoot + bHalf*tan(0) ,  -bHalf , zWing ];
    tipR.J = diag([8e-5 1.5e-4 2.2e-4]);
    tipL.J = tipR.J;
    
    % ======================================================================
    %  2. Combine
    % ======================================================================
    % lumps   = [fus bat tail fin tipR tipL];
    lumps   = [fus tail fin tipR tipL];
    lumpNoWing = [fus tail fin];
    mTot    = sum([lumps.m]);
    rCM     = sum( cat(1,lumps.m).*cat(1,lumps.r) , 1)/mTot;
    
    Ishift  = @(J,r,m) J + m*(dot(r,r)*eye(3) - r(:)*r(:).');     % -- parallel axis
    
    Jtot = zeros(3);
    for L = lumps
        Jtot = Jtot + Ishift(L.J , L.r - rCM , L.m);
    end
    
    % ======================================================================
    %  3. Output structure
    % ======================================================================
    params.mass  = mTot;
    params.mNoWing = sum([lumpNoWing.m]);
    params.rCM   = rCM(:);
    params.J     = Jtot;
    
    % keep sub-parts if caller wants them
    params.components = struct('fuselage',fus,'battery',bat,'tail',tail, ...
                            'tipR',tipR,'tipL',tipL);
    
    params.components.fuselage.M6comp_rCM = [ fus.m*eye(3) , -fus.m*skew(rCM);
               fus.m*skew(rCM) , fus.J ]; % Not sure if to use rCM or just r

    params.components.fuselage.M6comp_r = [ fus.m*eye(3) , -fus.m*skew(fus.r);
               fus.m*skew(fus.r) , fus.J ];
    % convenient 6×6 spatial mass
    params.M6 = [ mTot*eye(3) , -mTot*skew(rCM);
               mTot*skew(rCM) , Jtot ];
end
% 
% function prm = paramsRigid_PazyUAV()
% % paramsRigid_PazyUAV   Rigid–body lumped-parameter model for Pazy UAV
% %
% % All positions expressed in *body axes* (X forward, Y starboard, Z down)
% % Units: metres, kilograms, kg·m²
% % ----------------------------------------------------------------------
% 
% % ----- user-tunable geometric helpers ---------------------------------
% c      = 0.10;             % wing chord [m]
% ea     = 0.42*c;           % elastic-axis x-offset from LE
% bHalf  = 0.55;             % half-span to tip mass
% xWingRoot = 0.0;           % place root at X = 0
% zWing     = 0.0;           % mid-plane
% 
% % ======================================================================
% %  1. Lumps
% % ======================================================================
% 
% % Fuselage  -------------------------------------------------------------
% fus.m  = 0.60;
% fus.r  = [ xWingRoot-ea/2 ,   0 , zWing ];      % put CM slightly ahead of EA
% fus.J  = diag([7.1e-3  , 2.3e-2 , 2.3e-2]);
% 
% % Battery (inside fuse) -------------------------------------------------
% bat.m  = 0.18;
% bat.r  = fus.r + [0.08 , 0 , 0];                % 8 cm forward
% bat.J  = diag([3.3e-4 , 1.9e-3 , 2.2e-3]);
% 
% % Horizontal tail  ------------------------------------------------------
% tail.m = 0.12;
% tail.r = fus.r + [-0.60 , 0 , 0];               % 0.6 m aft
% tail.J = diag([2.0e-3 , 1.1e-3 , 3.0e-3]);
% 
% % Wing tip masses (actuators) ------------------------------------------
% tipR.m = 0.035;  tipL.m = tipR.m;
% tipR.r = [ xWingRoot + bHalf*tan(0) ,  +bHalf , zWing ];
% tipL.r = [ xWingRoot + bHalf*tan(0) ,  -bHalf , zWing ];
% tipR.J = diag([8e-5 1.5e-4 2.2e-4]);
% tipL.J = tipR.J;
% 
% % ======================================================================
% %  2. Combine
% % ======================================================================
% lumps   = [fus bat tail tipR tipL];
% mTot    = sum([lumps.m]);
% rCM     = sum( cat(1,lumps.m).*cat(1,lumps.r) , 1)/mTot;
% 
% Ishift  = @(J,r,m) J + m*(dot(r,r)*eye(3) - r(:)*r(:).');     % -- parallel axis
% 
% Jtot = zeros(3);
% for L = lumps
%     Jtot = Jtot + Ishift(L.J , L.r - rCM , L.m);
% end
% 
% % ======================================================================
% %  3. Output structure
% % ======================================================================
% prm.mass  = mTot;
% prm.rCM   = rCM(:);
% prm.J     = Jtot;
% 
% % keep sub-parts if caller wants them
% prm.components = struct('fuselage',fus,'battery',bat,'tail',tail, ...
%                         'tipR',tipR,'tipL',tipL);
% 
% % convenient 6×6 spatial mass
% prm.M6 = [ mTot*eye(3) , -mTot*skew(rCM);
%            mTot*skew(rCM) , Jtot ];
% 
% % Fuse cylinder sample (η from 0 at nose to 1 at tail)
% eta = linspace(0,1,30);
% R   = 0.04;                     % 4 cm radius
% Xc  = xWingRoot - ea/2  + 0.8*eta;      % 0.8 m long
% fuseCoords = [ Xc ;
%                zeros(size(Xc)) ;
%                R*cos(2*pi*eta) ];        % simple sinusoid
% prm.geom.fuselageLine = fuseCoords;
% end
% 
% function S = skew(v), S = [  0   -v(3)  v(2);
%                              v(3)   0   -v(1);
%                             -v(2) v(1)    0 ]; end
% 
% % function prm = paramsRigid_PazyUAV()
% % % Rigid-body parameters for the "Pazy-UAV" test case
% % % --------------------------------------------------
% % m_fuse   = 0.60;  
% % r_fuse = [-0.10  0   0];
% % Ifuse    = diag([7.1e-3 2.3e-2 2.3e-2]);
% % prm.fuse.m_fuse = m_fuse;prm.fuse.r_fuse = r_fuse;prm.fuse.Ifuse = Ifuse;
% % 
% % m_tail   = 0.12;  
% % % r_tail = [-0.70  0   0];
% % r_tail = [-0.70  0   0];
% % Itail    = diag([2.0e-3 1.1e-3 3.0e-3]);
% % 
% % m_batt   = 0.18;  
% % r_batt = [-0.10  0   0];
% % Ibatt    = diag([3.3e-4 1.9e-3 2.2e-3]);
% % 
% % m_wtip   = 0.035; 
% % r_wR   = [ 0.60  0.55 0];
% % I_wR     = diag([8.0e-5 1.5e-4 2.2e-4]);
% % 
% % m_wtipL  = m_wtip; 
% % r_wL   = [ 0.60 -0.55 0];
% % I_wL     = I_wR;
% % 
% % % ---- 1. Lump to total CM -------------------------------------------
% % mTot = m_fuse+m_tail+m_batt+2*m_wtip;
% % rCM  = (m_fuse*r_fuse + m_tail*r_tail + m_batt*r_batt + ...
% %         m_wtip*r_wR   + m_wtipL*r_wL)/mTot;
% % 
% % prm.mass = mTot;
% % prm.rCM  = rCM(:);
% % 
% % % ---- 2. Parallel-axis: shift each inertia to the total CM ----------
% % I_CM = @(I,r,m) I + m*( dot(r,r)*eye(3) - r(:)*r(:).' );
% % 
% % J   = I_CM(Ifuse ,r_fuse-rCM,m_fuse) ...
% %     + I_CM(Itail ,r_tail-rCM,m_tail) ...
% %     + I_CM(Ibatt ,r_batt-rCM,m_batt) ...
% %     + I_CM(I_wR  ,r_wR -rCM,m_wtip ) ...
% %     + I_CM(I_wL  ,r_wL -rCM,m_wtip );
% % 
% % prm.J = J;
% % 
% % % ---- 3. 6×6 spatial mass matrix ------------------------------------
% % prm.M6 = [ mTot*eye(3) , -mTot*skew(rCM);
% %            mTot*skew(rCM), J ];
% % 
% % end
% % 
function S = skew(v)
S=[ 0   -v(3)  v(2);
    v(3)    0  -v(1);
   -v(2) v(1)     0 ];
end
