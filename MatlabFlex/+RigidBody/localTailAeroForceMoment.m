
function [F_tail, M_tail, aux] = localTailAeroForceMoment(U, alpha, delta_e, cfg, q_pitch)
%LOCALTAILAEROFORCEMOMENT Horizontal-tail force and moment in body axes.
%
% Body-axis convention:
%   x_B forward
%   y_B right
%   z_B down
%
% Inputs:
%   U        : airspeed at the aircraft reference point [m/s]
%   alpha    : body angle of attack, atan2(w,u) [rad]
%   delta_e  : elevator deflection [rad]
%   cfg      : configuration structure
%   q_pitch  : optional body pitch rate q [rad/s]
%
% Outputs:
%   F_tail   : tail aerodynamic force in body axes [N]
%   M_tail   : tail aerodynamic moment about the rigid-body CG [N*m]
%   aux      : diagnostic structure
%
% Notes:
%   The moment is always computed as
%
%       M_tail = M_ac + cross(r_tail_from_CG, F_tail)
%
%   This is deliberately not replaced by +/- Lt*l_t shortcuts.  Those
%   shortcuts are easy to get wrong when switching between z-up and z-down
%   conventions.

if nargin < 5 || isempty(q_pitch)
    q_pitch = 0;
end

rho = localGetNested(cfg, {'flight','rho'}, 1.225);

S_t  = localGetNested(cfg, {'tail','S'}, ...
       localGetNested(cfg, {'tail','area'}, 0.035));
CL_a = localGetNested(cfg, {'tail','CL_alpha'}, 4.4);
eta  = localGetNested(cfg, {'tail','eta_delta'}, 1.0);

i_t = deg2rad(localGetNested(cfg, {'tail','incidence_deg'}, 0.01));

CD0 = localGetNested(cfg, {'tail','CD0'}, 0.01);
e   = localGetNested(cfg, {'tail','e'}, 0.8);
AR  = localGetNested(cfg, {'tail','AR'}, 4.0);

c_t = localGetNested(cfg, {'tail','c'}, ...
      localGetNested(cfg, {'tail','chord'}, 0.10));

Cm0 = localGetNested(cfg, {'tail','Cm0'}, 0.0);
Cm_delta = localGetNested(cfg, {'tail','Cm_delta'}, 0.0);

l_t = localGetNested(cfg, {'tail','arm'}, 0.60);

% Moment arm from CG to tail aerodynamic center in body axes.
% Prefer an explicit vector.  Fall back to an aft tail: negative x_B.
if isfield(cfg,'tail') && isfield(cfg.tail,'r_B') && ~isempty(cfg.tail.r_B)
    rTail_B = cfg.tail.r_B(:);
elseif isfield(cfg,'rigidEOMset') && isfield(cfg.rigidEOMset,'rTail_B') && ...
        ~isempty(cfg.rigidEOMset.rTail_B)
    rTail_B = cfg.rigidEOMset.rTail_B(:);
elseif isfield(cfg,'geom') && isfield(cfg.geom,'rTail_B') && ~isempty(cfg.geom.rTail_B)
    rTail_B = cfg.geom.rTail_B(:);
else
    rTail_B = [-l_t; 0; 0];
end

if numel(rTail_B) ~= 3
    error('localTailAeroForceMoment:BadTailArm', ...
        'cfg.tail.r_B must be a 3-vector from CG to tail aerodynamic center.');
end

Ueff = max(abs(U), 1e-9);
qbar = 0.5*rho*Ueff^2;

% Include the local velocity induced by pitch rate at the tail.
% For an aft tail r_x < 0 and positive q, omega x r gives positive w_B
% at the tail.  In z-down axes this increases local alpha.
vRef_B  = [Ueff*cos(alpha); 0; Ueff*sin(alpha)];
omega_B = [0; q_pitch; 0];
vTail_B = vRef_B + cross(omega_B, rTail_B);

alphaTailKin = atan2(vTail_B(3), vTail_B(1));

alphaEff = alphaTailKin - i_t - eta*delta_e;

CLt = CL_a*alphaEff;
CDt = CD0 + CLt^2/(pi*e*AR);

Lt = qbar*S_t*CLt;
Dt = qbar*S_t*CDt;

% Keep the same simple body-axis force convention used in the existing trim:
% drag opposes +x_B and lift acts upward, i.e. -z_B for positive CL.
F_tail = [-Dt; 0; -Lt];

M_ac_B = [0; qbar*S_t*c_t*(Cm0 + Cm_delta*delta_e); 0];

M_tail = M_ac_B + cross(rTail_B, F_tail);

aux = struct();
aux.qbar = qbar;
aux.rTail_B = rTail_B;
aux.alphaTailKin = alphaTailKin;
aux.alphaEff = alphaEff;
aux.CL = CLt;
aux.CD = CDt;
aux.L = Lt;
aux.D = Dt;
aux.M_ac_B = M_ac_B;
end

function val = localGetNested(S, fields, defaultVal)
val = defaultVal;
try
    tmp = S;
    for i = 1:numel(fields)
        if ~isstruct(tmp) || ~isfield(tmp, fields{i}) || isempty(tmp.(fields{i}))
            return
        end
        tmp = tmp.(fields{i});
    end
    val = tmp;
catch
    val = defaultVal;
end
end


% function [F_tail, M_tail] = localTailAeroForceMoment(U, alpha, delta_e, cfg)
% 
% rho = localGet(cfg.flight,'rho',1.225);
% q = 0.5*rho*U^2;
% 
% S_t = localGetNested(cfg, {'tail','S'}, localGetNested(cfg, {'tail','area'}, 0.035));
% CL_a = localGetNested(cfg, {'tail','CL_alpha'}, 4.4);
% eta = localGetNested(cfg, {'tail','eta_delta'}, 1.0);
% i_t = deg2rad(localGetNested(cfg, {'tail','incidence_deg'}, 0.01));
% CD0 = localGetNested(cfg, {'tail','CD0'}, 0.01);
% e = localGetNested(cfg, {'tail','e'}, 0.8);
% AR = localGetNested(cfg, {'tail','AR'}, 4.0);
% l_t = localGetNested(cfg, {'tail','arm'}, 0.6);
% 
% CLt = CL_a*(alpha - i_t - eta*delta_e);
% CDt = CD0 + CLt^2/(pi*e*AR);
% 
% Lt = q*S_t*CLt;
% Dt = q*S_t*CDt;
% 
% F_tail = [-Dt; 0; -Lt];
% 
% if isfield(cfg,'tail') && isfield(cfg.tail,'r_B')
%     rTail = cfg.tail.r_B(:);
%     M_tail = cross(rTail,F_tail);
% else
%     % Preserve the sign convention used in the existing routine.
%     M_tail = [0; Lt*l_t; 0];
% end
% end
% 
% function val = localGet(S, field, defaultVal)
% if isstruct(S) && isfield(S,field) && ~isempty(S.(field))
%     val = S.(field);
% else
%     val = defaultVal;
% end
% end
% function val = localGetNested(S, fields, defaultVal)
% val = defaultVal;
% try
%     tmp = S;
%     for i = 1:numel(fields)
%         if ~isstruct(tmp) || ~isfield(tmp, fields{i}) || isempty(tmp.(fields{i}))
%             return
%         end
%         tmp = tmp.(fields{i});
%     end
%     val = tmp;
% catch
%     val = defaultVal;
% end
% end