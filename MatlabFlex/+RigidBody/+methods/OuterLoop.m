classdef OuterLoop < AeroFlex.ctrl.ControllerBase
    % Thin wrapper that sits above the wing–level NMPC/MHE
    properties
        eom     RigidBody.EOM
        ctrl    RigidBody.LQRController
        Ts
    end
    methods
        function obj = OuterLoop(cfg,wingTS)
            obj@AeroFlex.ctrl.ControllerBase(cfg);
            obj.Ts  = cfg.ctrl.TsOuter;
            obj.eom = RigidBody.EOM(cfg.uav);
            % --- small-signal linearisation (around trim straight flight)
            x0  = StateQuat.eyeQuat;             % attitude
            x0  = [x0(:); zeros(9,1)];           % all zeros else
            A   = zeros(6);  B = eye(6);         %# placeholder (supply real)
            Q   = diag([2 2 1   1 1 1]);         % angles & vels
            R   = eye(6)*0.05;
            obj.ctrl = RigidBody.LQRController(A,B,Q,R);
        end
        %---------------------------------------------------------------
        function tau = computeControl(obj,xrbSmall)
            tau = obj.ctrl.compute(xrbSmall);
        end
    end
end
