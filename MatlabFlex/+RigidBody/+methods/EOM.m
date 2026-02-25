% classdef EOM < handle
%     % Quaternion-based 6-DOF rigid-body equations in body axes.
%     % State: x = [ q0 q1 q2 q3  u v w  p q r  x y z ]'  (4 + 3 + 3 + 3)
%     properties
%         mass   double
%         Inertia double    % 3x3 body inertia
%         InvI   double
%         gVec   double = [0;0;9.807];   % +Z down (NED)
%     end
%     methods
%         function obj = EOM(cfg, prms)
%             obj.mass    = prms.mass;
%             obj.Inertia = prms.J;
%             obj.InvI    = inv(obj.Inertia);
%         end
% 
%         function xdot = f(obj, x, uB, Fb, Mb)
%             % Inputs are in body frame.
%             q  = x(1:4);      Vb = x(5:7);   Om = x(8:10);   rI = x(11:13); 
%             Rb2I = quat2dcm(q.');                   % body->inertial
%             Fg_b = Rb2I.'*(obj.mass*obj.gVec);      % weight in body axes (+Z down)
% 
%             acc_b  = (uB + Fb + Fg_b)/obj.mass - cross(Om, Vb);
%             Om_dot = obj.InvI*( Mb - cross(Om, obj.Inertia*Om) );
% 
%             Omega  = [ 0, -Om.' ; Om, -skew(Om) ];
%             q_dot  = 0.5*Omega*q;
% 
%             rI_dot = Rb2I.'*Vb;
% 
%             xdot = [q_dot; acc_b; Om_dot; rI_dot];
%         end
% 
%         function [xNew,STM] = step(obj,x,uExt,dt)
%             % 4th-order RK with quaternion renormalization
%             F_b = uExt.F;   M_b = uExt.M;   uB = uExt.uB;
%             k1 = obj.f(x, uB, F_b, M_b);
%             k2 = obj.f(x + 0.5*dt*k1, uB, F_b, M_b);
%             k3 = obj.f(x + 0.5*dt*k2, uB, F_b, M_b);
%             k4 = obj.f(x + dt*k3,     uB, F_b, M_b);
%             xNew = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);
% 
%             % Normalize quaternion (critical for stability)
%             xNew(1:4) = xNew(1:4)/norm(xNew(1:4));
%             STM = [];   % placeholder
%         end
%     end
% end
% 
% function S = skew(v)
% S = [  0   -v(3)  v(2);
%       v(3)   0   -v(1);
%      -v(2)  v(1)   0 ];
% end

classdef EOM < handle
    % Quaternion-based 6-DOF rigid-body equations *in body frame*
    % x = [ q0 q1 q2 q3  u v w  p q r ]'    (4 + 3 + 3)
    %        quat      vel     ang-vel
    % External forces:  F_b  (N) and M_b (N·m) passed from FlexCoupler
    % 6-DOF rigid-body equations (quaternion form)
    %   Uses body axes, wind‐axes aerodynamic forces come from AeroFlex ROM
    %
    properties
        mass   double
        Inertia double  % 3×3 matrix in body frame
        InvI   double
        gVec   double = [0 0 9.807].';   % +Z down (NED)
    end
    methods
        function obj = EOM(cfg, prms)
            obj.mass   = prms.mass;
            % obj.Inertia= prms.M6;
            obj.Inertia= prms.J;
            cfg.Ibody = prms.J;
            % obj.InvI   = inv(cfg.Ibody);
            obj.InvI   = inv(cfg.Ibody);
        end
        %---------------------------------------------------------------
        function xdot = f(obj,t, x,uB,Fb,Mb)
            % x .. 13×1 StateQuat
            % uB .. thrust force in body frame  (3×1)
            % Fb .. aerodynamic forces in body frame (3×1)
            % Mb .. aerodynamic moments in body frame (3×1)
            q  = x(1:4);     Vb = x(5:7);   Om = x(8:10);   Pi = x(11:13);
            R  = quat2dcm(q.');             % body→inertial
            % translational accelerations
            Fi = R.'*obj.mass*obj.gVec;      % gravity in body frame
            acc_b = (uB + Fb + Fi)/obj.mass - cross(Om,Vb);
            % rotational accelerations
            Om_dot = obj.InvI * ( Mb - cross(Om,obj.Inertia*Om) );
            % quaternion kinematics
            Omega = [   0   -Om.';
                         Om   -skew(Om)];
            q_dot  = 0.5*Omega*q;
            % inertial position rate
            Pi_dot = R.'*Vb;
            % stacked
            xdot = [q_dot; acc_b; Om_dot; Pi_dot];
        end
        function [xNew,STM] = step(obj,x,uExt,dt)
            % RK4 with STM (optional)
            F_b = uExt.F;   M_b = uExt.M;
            uB = uExt.uB ;
            k1  = obj.f(0,x,uB, F_b,M_b);
            k2  = obj.f(0,x+0.5*dt*k1,uB, F_b,M_b);
            k3  = obj.f(0,x+0.5*dt*k2,uB, F_b,M_b);
            k4  = obj.f(0,x+dt*k3 ,uB, F_b,M_b);
            xNew= x + dt/6*(k1+2*k2+2*k3+k4);
            STM = eye(numel(x));     % Temp, neeed to implement
        end
    end
end

function S = skew(v)
S = [   0   -v(3)  v(2);
       v(3)   0   -v(1);
      -v(2)  v(1)   0 ];
end