classdef FlexCoupler < handle
    % Bridges AeroFlex wing ROM with rigid-body EoM
    properties
        wing   AeroFlex.sim.ROMIntegrator
        RB     RigidBody.EOM6DOFQuat
        LUT    struct
        dt double
    end
    methods
        function obj = FlexCoupler(cfg,rbMass,rbI)
            obj.wing = cfg.modelHandle(cfg);
            obj.RB   = RigidBody.EOM(rbMass,rbI);
            obj.dt   = cfg.sim.dt;
            obj.LUT  = load(fullfile('data','aeroGrid.mat')).aeroLUT;
        end
        function [xRB,xWing] = step(obj,xRB,xWing,uCmd,wgust)
            % 1) apply flap to wing ROM, integrate internal flexible modes
            [xWing,forcesWing] = obj.wing.step(xWing,uCmd.deltaDeg,wgust,[],false);
            % 2) map the aerodynamic coefficients to body-frame forces
            coeff = obj.interpAero(xWing,uCmd,xRB);
            F_b   = coeff.CL;   % expand to 3-axes!
            M_b   = coeff.CM;
            % 3) rigid body propagation
            [xRB,~] = obj.RB.step(xRB,struct('F',F_b,'M',M_b),obj.dt);
        end
    
    function coeff = interpAero(obj,xWing,uCmd,xRB)
        % Very light: tri-linear interp in α-δ-U -> CL,CD,CM
        alpha = xRB2alpha(xRB);   deltaDeg = uCmd.deltaDeg;  U = norm(xRB(5:7));
        coeff.CL = interpn(obj.LUT.alpha,obj.LUT.delta,obj.LUT.U,obj.LUT.CL,...
                           alpha,deltaDeg,U,'linear');
        coeff.CD = interpn(obj.LUT.alpha,obj.LUT.delta,obj.LUT.U,obj.LUT.CD,...
                           alpha,deltaDeg,U,'linear');
        coeff.CM = interpn(obj.LUT.alpha,obj.LUT.delta,obj.LUT.U,obj.LUT.CM,...
                           alpha,deltaDeg,U,'linear');
    end
    function a = xRB2alpha(xRB), a = atan2d(-xRB(6),xRB(5)); end
    end
end
