classdef ROMIntegrator < handle
    properties
        L, Lfac, Ufac, piv, gamma, delta, dt, idx, parConst
    end
    methods
        function obj = ROMIntegrator(cfg,beam,aero,base)
            obj.dt = cfg.sim.dt;                           % 1. matrices
            [obj.L,~] = AeroFlex.core.assemble_L_matrix(cfg,beam,aero,base);
            obj.gamma = (2-sqrt(2))/2;   obj.delta = -2*sqrt(2)/3;

            A = speye(size(obj.L)) - obj.gamma*obj.dt*obj.L;
            [obj.Lfac,obj.Ufac,obj.piv] = lu(A,'vector');

            obj.idx      = AeroFlex.core.buildIndexStruct(beam.Nm,aero.Na);
            % obj.parConst = buildParConst(beam,aero,base,cfg);  
            % forces0 = beam.eta_e + aero.DataMatrix.forces_aero_beam_dof;
            forces0 = beam.eta_e + aero.DataMatrix.forces_aero_beam_dof;
            % disp(beam.eta_e);
            % forces0 = aero.DataMatrix.forces_aero_beam_dof;
            % forces0 = beam.eta_e;
             %── build constant parameter sub‑struct -----------------------
            obj.parConst = struct( ...
                'Gamma1',   beam.Gamma1, ...
                'Gamma2',   beam.Gamma2, ...
                'Gamma_g',  base.Gamma_g, ...
                'Gamma_xi', base.Gamma_xi, ...
                'forces_0', forces0, ...   % steady aero load Nm×1
                'scaleAero', cfg.flight.a*cfg.flight.t_inf, ...
                't_inf', cfg.flight.t_inf, ...
                'Fscale', cfg.flight.Fscale, ...
                'scaleA', cfg.flight.a, ...
                'gustSet', cfg.gust, ...
                'Bw', aero.forceMap.Bw, ...
                'Dw', aero.forceMap.Dw, ...
                'Bdel',  aero.forceMap.B_delta, ...
                'Ddel',  aero.forceMap.D_delta, ...
                'Bddel', aero.forceMap.B_ddelta, ...
                'Dddel', aero.forceMap.D_ddelta, ...
                'gust_input', aero.gust_input, ...
                'dt', cfg.sim.dt, ...
                'Na', aero.Na, ...
                'SwTest', cfg.SwTest, ...
                'N_Thrust', cfg.trim.thrust);
        end
        %------------------------------------------------------------------
        function [x_np1,S_np1] = step(obj,x_n,u_k,g_k,S_n,storeSens)
            if nargin<6, storeSens=false;  end
            pc = obj.parConst;  pc.gust=g_k; 
            % pc.u=u_k;         % 2. params
            pc.u_ctrl=u_k;         % 2. params

            dt = obj.dt;   
            g  = obj.gamma;  
            d = obj.delta;     % stage 1
            
            k1N = AeroFlex.sim.nonlinear_terms(x_n,pc,obj.idx)*dt;
            if storeSens
                [Nq1,Nu1] = AeroFlex.sim.nonlinearJacobian(x_n,obj.idx,pc);
                j1N = dt*(Nq1*S_n + Nu1);
                rhsS= obj.L*(S_n+g*j1N)*dt;
                j1L = obj.solveImplicit(rhsS);
                S1  = S_n + g*(j1L+j1N);
            end
            k1L = obj.solveImplicit(obj.L*(x_n+g*k1N)*dt);
            y1  = x_n + g*(k1L+k1N);

            % stage 2
            k2N = AeroFlex.sim.nonlinear_terms(y1,pc,obj.idx)*dt;
            if storeSens
                [Nq2,Nu2] = AeroFlex.sim.nonlinearJacobian(y1,obj.idx,pc);
                j2N = dt*(Nq2*S1+Nu2);
                rhsS= obj.L*(S_n+(1-g)*j1L+d*j1N+(1-d)*j2N)*dt;
                j2L = obj.solveImplicit(rhsS);
                S2  = S_n+(1-g)*j1L+g*j2L+d*j1N+(1-d)*j2N;
            end
            k2L = obj.solveImplicit(obj.L*(x_n+(1-g)*k1L+d*k1N+(1-d)*k2N)*dt);

            % stage 3
            y2  = x_n+(1-g)*k1L+g*k2L+d*k1N+(1-d)*k2N;
            k3N = AeroFlex.sim.nonlinear_terms(y2,pc,obj.idx)*dt;
            if storeSens
                [Nq3,Nu3] = AeroFlex.sim.nonlinearJacobian(y2,obj.idx,pc);
                j3N = dt*(Nq3*S2+Nu3);
                S_np1 = S_n+(1-g)*(j1L+j2N)+g*(j2L+j3N);
            else
                S_np1 = [];
            end
            x_np1 = x_n + (1-g)*(k1L+k2N) + g*(k2L+k3N);
        end
        %------------------------------------------------------------------
        function y = solveImplicit(obj,rhs)
            z = rhs(obj.piv,:);   z = obj.Lfac \ z;   y = obj.Ufac \ z;
        end
    end
end
