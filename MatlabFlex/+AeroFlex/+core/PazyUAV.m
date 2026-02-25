classdef PazyUAV < handle
    %PazyUAV Class for visualizing the Pazy wing with attached rigid UAV components.
    % This class renders a 3D view of the flexible Pazy wing (using its FEM nodes 
    % and aerodynamic data) together with the rigid fuselage, tail, battery, and tip masses.
    % Usage:
    %   viz = PazyUAV(wingNodes, af);
    %   viz.plotStructure();
    %
    % Inputs:
    %   wingNodes : N x 3 array of [x,y,z] coordinates of the wing elastic axis nodes (e.g. quarter-chord line).
    %   af        : 2 x M array of airfoil coordinates [x/c; z/c] defining the airfoil shape (NACA0018) in normalized units.
    %              The airfoil should be defined from leading edge around to trailing edge (closed loop).
    %
    % The plotStructure method will create a new figure with the wing, fuselage (as a cylinder),
    % and markers for fuselage mass, tail mass, battery, and tip masses. The elastic axis, chord lines, 
    % and airfoil cross-section patches are drawn for clarity.
    
    properties
        wingNodes   % Nx3 array of wing elastic axis node coordinates (global coordinate system)
        airfoil    % 2xM array of airfoil shape (x/c, z/c)
        paramsRB   % structure of rigid-body parameters (masses, inertias, positions, fuselage geometry)
        chord      % (scalar) wing chord length (meters)
    end
    
    methods
        function obj = PazyUAV(wingNodes, af, chord)
            % Constructor for PazyUAV
            if nargin < 3
                chord = 0.10;  % Default chord length 0.10 m (100 mm) for Pazy wing
            end
            obj.wingNodes = wingNodes;
            obj.airfoil = af;
            obj.chord = chord;
            % Load rigid-body parameters (masses, CG locations, etc.)
            obj.paramsRB = paramsRigid_PazyUAV();
        end
        
        function plotStructure(obj)
            % plotStructure Render the 3D geometry of the flexible wing and rigid-body components.
            
            % Prepare figure and axes
            figure; 
            ax = gca;
            hold(ax, 'on');
            axis equal;  % equal scaling for proper geometry
            xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
            title('Pazy UAV Flexible Wing and Rigid-Body Visualization');
            
            % Extract wing node coordinates for convenience
            nodes = obj.wingNodes;
            % Ensure we have a column for each of x, y, z
            if size(nodes,2) ~= 3
                error('wingNodes must be an N x 3 array of [x,y,z] coordinates.');
            end
            
            % Determine the center (wing root) node index (assuming one at or near [x≈0,y=0,z≈0]):
            [~, centerInd] = min(sum(nodes.^2, 2));  % this finds the node closest to origin
            % Alternatively, find node with y ~ 0 (if center is at y=0):
            centerCandidates = find(abs(nodes(:,2)) < 1e-6);
            if ~isempty(centerCandidates)
                % if multiple near y=0, pick the one with smallest |x|+|z|
                [~, idx] = min(abs(nodes(centerCandidates,1)) + abs(nodes(centerCandidates,3)));
                centerInd = centerCandidates(idx);
            end
            
            % Separate left and right side node indices relative to center
            N = size(nodes,1);
            leftIndices = 1:centerInd;   % include center in left side
            rightIndices = centerInd:N;  % include center in right side
            % Ensure center is not plotted twice for certain elements
            % We'll handle center separately when plotting cross-sections to avoid duplication
            
            % Plot elastic axis (quarter-chord line) through all wing nodes
            plot3(nodes(:,1), nodes(:,2), nodes(:,3), 'k--', 'LineWidth', 1.5);  % dashed black line for elastic axis
            
            % Precompute the airfoil shape coordinates scaled to chord, with quarter-chord at origin
            af = obj.airfoil;
            % Shift airfoil so that quarter-chord (x/c = 0.25) is at (0,0)
            af_coords = af;
            af_coords(1,:) = af(1,:) - 0.25;  % now leading edge ~ +0.25, trailing edge ~ 0.75 from origin
            % Scale to actual chord length
            af_coords = af_coords * obj.chord;
            
            % Decide on which spanwise stations to draw chord lines and airfoil patches
            stations = [];  % indices of wingNodes to draw cross-sections at
            if ~isempty(centerInd)
                stations(end+1) = centerInd;
            end
            % Mid-span (approximately halfway between root and tip) for each side if enough nodes
            if centerInd > 2
                % Left side mid (exclude center and tip)
                midLeftIndex = floor((1 + (centerInd)) / 2);
                if midLeftIndex < centerInd && midLeftIndex > 1
                    stations(end+1) = midLeftIndex;
                end
            end
            if N - centerInd > 1
                % Right side mid
                midRightIndex = floor((centerInd + 1 + N) / 2);
                if midRightIndex > centerInd && midRightIndex < N
                    stations(end+1) = midRightIndex;
                end
            end
            % Wing tips (both ends)
            stations(end+1) = 1;    % left tip (assuming index1 is left tip)
            stations(end+1) = N;    % right tip
            
            stations = unique(stations);  % remove any duplicates and sort
            % Plot chord lines and airfoil cross-section patches at the selected stations
            for ii = 1:length(stations)
                idx = stations(ii);
                pos = nodes(idx, :);  % quarter-chord position at this station
               % For straight wing, we can directly use global axes for X and Z.
                Xcoords = pos(1) + af_coords(1,:);   % add quarter-chord X position
                Zcoords = pos(3) + af_coords(2,:);   % add vertical position (Z)
                % Y coordinates: all the airfoil points have the same spanwise coordinate as the node
                Ycoords = pos(2) * ones(size(Xcoords));
                
                % Draw the chord line (from leading edge to trailing edge) at this station
                % Leading edge is quarter-chord + 0.25*c (forward), trailing edge is quarter-chord - 0.75*c (aft)
                LE_global = [pos(1) + 0.25*obj.chord, pos(2), pos(3)];
                TE_global = [pos(1) - 0.75*obj.chord, pos(2), pos(3)];
                plot3([LE_global(1), TE_global(1)], [LE_global(2), TE_global(2)], [LE_global(3), TE_global(3)], ...
                      'b-', 'LineWidth', 1.2);
                
                % Draw the airfoil cross-section shape as a filled patch
                patch('XData', Xcoords, 'YData', Ycoords, 'ZData', Zcoords, ...
                      'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k', 'FaceAlpha', 1.0);
            end
            
            %% Plot the fuselage as a cylinder centered at the wing root
            fus_R = obj.paramsRB.fuselageGeom.radius;
            fus_L = obj.paramsRB.fuselageGeom.length;
            % Create a cylinder aligned with X-axis (nose-to-tail)
            [Xc, Yc, Zc] = cylinder(fus_R, 20);  % default cylinder along Z from 0 to 1
            % Transform to align along X-axis from -L/2 to +L/2 (centered at 0)
            Xcoords = Zc * fus_L - fus_L/2;   % scale and shift to length L, centered at 0
            Ycoords = Xc;                    % original X circle -> Y cross-section
            Zcoords = Yc;                    % original Y circle -> Z cross-section
            % Plot cylinder surface
            surf(ax, Xcoords, Ycoords, Zcoords, 'FaceColor', [0.4 0.6 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            
            % Plot mass components as small spheres (or markers) at specified locations
            % Choose a small radius for visualization of masses
            sphereRadius = 0.01;  % 1 cm radius for mass markers
            [sx, sy, sz] = sphere(10);  % sphere mesh for marker
            sx = sphereRadius * sx; sy = sphereRadius * sy; sz = sphereRadius * sz;
            
            % Fuselage mass (assumed at fuselage CG position)
            fusCG = obj.paramsRB.fuselage.pos;  % should be [x,y,z] relative to wing root
            h1 = surf(ax, sx + fusCG(1), sy + fusCG(2), sz + fusCG(3));
            set(h1, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.9);  % black sphere
            
            % Battery mass
            batPos = obj.paramsRB.battery.pos;
            h2 = surf(ax, sx + batPos(1), sy + batPos(2), sz + batPos(3));
            set(h2, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.9);  % red sphere
            
            % Tail mass
            tailPos = obj.paramsRB.tail.pos;
            h3 = surf(ax, sx + tailPos(1), sy + tailPos(2), sz + tailPos(3));
            set(h3, 'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.9);  % green sphere
            
            % Tip masses (left and right)
            tipPos = obj.paramsRB.tip.pos;  % should be 2x3 array: [pos_left; pos_right]
            if size(tipPos,1) == 2
                % Left tip mass
                h4 = surf(ax, sx + tipPos(1,1), sy + tipPos(1,2), sz + tipPos(1,3));
                set(h4, 'FaceColor', [1 0 1], 'EdgeColor', 'none', 'FaceAlpha', 0.9);  % magenta sphere
                % Right tip mass
                h5 = surf(ax, sx + tipPos(2,1), sy + tipPos(2,2), sz + tipPos(2,3));
                set(h5, 'FaceColor', [1 0 1], 'EdgeColor', 'none', 'FaceAlpha', 0.9);
            else
                % If only one set of tip coordinates provided, assume symmetric and mirror in Y
                posTip = tipPos(1,:);
                % left
                h4 = surf(ax, sx + posTip(1), sy + (-abs(posTip(2))), sz + posTip(3));
                set(h4, 'FaceColor', [1 0 1], 'EdgeColor', 'none', 'FaceAlpha', 0.9);
                % right
                h5 = surf(ax, sx + posTip(1), sy + abs(posTip(2)), sz + posTip(3));
                set(h5, 'FaceColor', [1 0 1], 'EdgeColor', 'none', 'FaceAlpha', 0.9);
            end
            
            % Finalize axes
            view(30, 15);  % set a good 3D view angle (azimuth, elevation)
            grid on;
            hold(ax, 'off');
        end
    end
end