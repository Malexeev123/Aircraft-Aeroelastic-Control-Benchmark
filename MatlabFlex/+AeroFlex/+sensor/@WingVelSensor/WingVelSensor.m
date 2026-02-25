% sensors/@WingVelSensor/WingVelSensor.m
classdef WingVelSensor < handle
% Returns v_y at five span‑wise stations 0,0.25,0.5,0.75,1 · L
% Compatible with modal state vector q1 (velocities).

    properties
        PhiY          % modal shape in Y direction (Ny×Nm)
        stations      % ξ/L positions (1×5)
    end

    methods
        function obj = WingVelSensor(beam, cfg)
            interp_node = []; L_interp = []; selIterp = [];
            % obj.stations = [0  0.25 0.5 0.75  1];
            obj.stations = [0.25 0.5 0.75  1];
            % obj.stations = [0.2 0.4 0.6 0.8  1];
            L_select = obj.stations*cfg.struct.L-cfg.struct.L/2;
            phi1_meas = zeros(6*length(L_select), beam.Nm);

            for i = 1:length(L_select)
                [I, ~] = find(beam.fem.coordinates(:,2) == L_select(i));
                selIdx = double((i-1)*6)+(1:6);
                if isempty(I)
                    selIterp(end+1,:) = selIdx;
                    L_interp(end+1) = L_select(i);

                else
                    phiIdx = double((I-1)*6)+(1:6);
                    phi_nodes = AeroFlex.core.localDofIndex(beam.fem.keepDofs, phiIdx);
                    if isempty(phi_nodes)
                        % interp_node(end+1) = I; 
                        L_interp(end+1) = L_select(i);
                        selIterp(end+1,:) = selIdx;
                       continue; 
                    end

                    phi1_meas(selIdx, :) = beam.phi1(phi_nodes, :);

                end
              
            end
            phi_interp = obj.interpModalShapes(beam,'y',L_interp, selIterp);
            obj.PhiY = phi1_meas+phi_interp;
        end

        function y = measure(obj,q1)  % Need to add the quaternion rotation
            y = obj.PhiY * q1;                   % v = Φ_y q̇
        end

        function H = jacobian(obj,~)  % need to redo this one
            H = obj.PhiY;    % linear in velocities
        end


        function phi_interp = interpModalShapes(obj, beam,rowIdx, s_interp, s_node)
            coords = beam.fem.coordinates(:,2);
            phi_interp = zeros(6*length(obj.stations), beam.Nm);
            for j = 1:length(s_interp)
                sA = s_interp(j);
                snodes = s_node(j,:);
                % --- locate the element that contains s --------------------------
                if abs(sA-coords(1)) < 1e-4               % exactly at the clamp
                    sL = coords(2);
                    sR = coords(end);
                    phi1L = beam.phi1(1:6, :);
                    phi1R = beam.phi1(end-5:end, :);
                else
                    iL = find(coords <= sA,1,'last');    % left node
                    if iL == numel(coords)              % s at wing tip -> clamp to tip
                        iL = iL-1;
                    end
                    iR = iL + 1;                        % right node
                    sL = coords(iL);
                    sR = coords(iR);
                    % xi = (sA - coords(sL)) / (coords(sR) - coords(sL));  % barycentric
                    % map = @(node) double((node-2)*6) + rowIdx(:);   % column vector of row indices
                    map = @(node) double((node-2)*6) + (1:6);   % column vector of row indices

                    rowsL = map(iL);    rowsR = map(iR);
                    phi1L  = beam.phi1(rowsL ,:);                % length(rowIdx)×Nm
                    phi1R  = beam.phi1(rowsR ,:);
                end
                elemLen = abs(sR - sL);
                phi_interp(snodes,:) = (phi1L*(sR - sA)/elemLen+ phi1R*(sA-sL)/elemLen);
            end
        end
    end
end
