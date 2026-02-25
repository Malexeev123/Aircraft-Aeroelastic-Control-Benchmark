classdef plotDOFSandActiveNodes
    %PLOTDOFSANDACTIVENODES Summary of this class goes here    
    properties
        fem
        rgParam
        paramsRB   
        aeromesh
        colWing  = [0.05 0.35 0.85]  % RGB for wing surface
    end
    
    methods
        function obj = plotDOFSandActiveNodes(beam, aero, rgParam)
            %PLOTDOFSANDACTIVENODES Construct an instance of this class
            obj.fem = beam.fem;
            obj.rgParam = rgParam;
            obj.aeromesh = aero.aeroMesh;
            % obj.plotNodes3D();
            % obj.plotNodes3D_V2();
            obj.paramsRB = RigidBody.methods.paramsRigid_PazyUAV();


        end
  
        function plotNodes3D_V2(obj)
            fem = obj.fem;
            figure; hold on; axis equal;
            xlabel('X'); ylabel('Y'); zlabel('Z');
            title('3D Node Layout');
        
            coords = fem.coordinates; % typically [num_node x 3] or [3 x num_node]
            if size(coords,2)==3 && size(coords,1)==fem.num_node
                % okay [node_i, (x,y,z)]
            elseif size(coords,1)==3 && size(coords,2)==fem.num_node
                coords = coords';
            end
            
            plot3(coords(:,1), coords(:,2), coords(:,3), 'ro','MarkerFaceColor','r');
            for i=1:fem.num_node
                text(coords(i,1), coords(i,2), coords(i,3), sprintf('Node %d', i));
            end
            obj.plotElements3D();
             hold off;
        end
        
        function plotElements3D(obj)
            % figure; hold on; axis equal;
            fem = obj.fem;
            coords = fem.coordinates;
            if size(coords,2)==3, coords=coords'; end
            
            for e=1:fem.num_elem
                nodesE = fem.connectivities(e,:);  % e.g. [n1,n2,n3]
                xyz = coords(:,nodesE);
                % For 3-node beams, you might just plot lines from node1->node2->node3:
                plot3(xyz(1,[1,3,2]), xyz(2,[1,3,2]), xyz(3,[1,3,2]), '-b','LineWidth',1.5);
            end
        end
        
        
        function plotNodes3D(obj)
        % plot the node coordinates with boundary colors
        fem = obj.fem;

        figure('Name','Pazy Nodes 3D','Color','white'); hold on; axis equal
        XYZ = fem.coordinates;
        bcs = fem.boundary_conditions;
        for n=1:fem.num_node
           x= XYZ(n,1); y=XYZ(n,2); z=XYZ(n,3);
           switch bcs(n)
             case 1, clr='r';  % clamp
             case -1,clr='g';  % tip
             otherwise, clr='b';
           end
           plot3(x,y,z,'o','MarkerFaceColor',clr,'MarkerSize',5,'MarkerEdgeColor','k');
           text(x,y,z, sprintf('%d',n),'Color','k','FontSize',9);
           
        end
        % also draw lines for each element
        ne = fem.num_elem;
        for e=1:ne
           c= fem.connectivities(e,:);
           cXYZ= XYZ(c,:);
           dx = cXYZ(:,1);
           dy = cXYZ(:,2);
           dz = cXYZ(:,3);
           plot3(cXYZ([1,3, 2],1), cXYZ([1,3, 2],2), cXYZ([1,3, 2],3), '-k','LineWidth',1.2);
        end
        xlabel('X'); ylabel('Y'); zlabel('Z');
        title('Nodes, boundary conditions, and connectivity');
        view(3);
        grid on;
        hold off;
        end

        function plotPazySystem(obj)
            %PLOTPAZYSYSTEM  Draw Pazy flexible wing mesh, idealised chord, and rigid
            % body masses (fuselage, tail, batteries, tip masses).
            %
            % Inputs
            %   beam   – instance of AeroFlex.model.Beam (or beam.fem struct)
            %   rgParam– struct from RigidBody.methods.paramsRigid_PazyUAV()
            %   chord  – physical chord length [m]   (default 0.10)
            %
            % Example
            %   cfg   = nominalConfig;
            %   beam  = cfg.modelHandle(cfg,beam,aero,base);  % your usual call
            %   rgPrm = RigidBody.methods.paramsRigid_PazyUAV();
            %   AeroFlex.vis.plotPazySystem(beam,rgPrm,0.10);
            
            % if nargin<3, chord = 0.10; end
            chord = obj.aeromesh.chord(1,1);
            fem = obj.fem;
            rgParam = obj.rgParam;
            XYZ = fem.coordinates;                    % [numNode × 3]
            if size(XYZ,2)~=3, XYZ = XYZ.'; end       % ensure N×3
            
            figure('Color','w','Name','Pazy-wing + Rigid Body'); clf; hold on; grid on
            axis equal; view(35,20);
            xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
            title('Pazy flexible wing + fuselage/tail point masses');
            
            %% --- 1. BEAM centre-line & nodes -----------------------------------
            plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'b.-','LineWidth',1.2,'MarkerSize',10);
            text(XYZ(:,1),XYZ(:,2),XYZ(:,3), string(1:size(XYZ,1)),'FontSize',8);
            
            %% --- 2. Draw a simple rectangular planform strip for each element ---
            % ea = 0.42*chord;                      % elastic axis ~42 % c
            ea = obj.aeromesh.elastic_axis(1,1);
            for e = 1:fem.num_elem
                nIdx  = fem.connectivities(e,[1 3 2]);      % [root mid tip] order
                pts   = XYZ(nIdx,:);                        % 3×3
                % Span-wise vector
                sVec  = pts(3,:)-pts(1,:);
                % Local chord direction (x-body forward) – here pick +x
                cVec  = [chord 0 0];
                % Leading & trailing-edge quad (flat plate)
                le    = pts + (ea)*cVec;   % leading edge
                te    = pts - (1-ea)*cVec; % trailing edge
                quadX = [le(:,1); flipud(te(:,1))];
                quadY = [le(:,2); flipud(te(:,2))];
                quadZ = [le(:,3); flipud(te(:,3))];
                patch(quadX,quadZ,quadY, 0.8*[1 1 1], 'FaceAlpha',0.3,'EdgeColor','none');
                % plot(quadX,quadY,quadZ, 0.8*[1 1 1], 'FaceAlpha',0.3,'EdgeColor','none');
            end
            
            %% --- 3. Plot fuselage / tail lumped-masses -------------------------
            plotMass(rgParam.rCM,  rgParam.mass,   'k','CM');
            addMass(rgParam, 'r_fuse', m_fuse,  'Fuse');
            addMass(rgParam, 'r_tail', m_tail,  'Tail');
            addMass(rgParam, 'r_batt', m_batt,  'Batt');
            addMass(rgParam, 'r_wR',   m_wtip,  'Tip-R');
            addMass(rgParam, 'r_wL',   m_wtipL, 'Tip-L');
            
            hold off
            
            %% nested helpers -----------------------------------------------------
                function plotMass(r, m, clr, lbl)
                    scatter3(r(1),r(2),r(3), 60,clr,'filled');
                    text(r(1),r(2),r(3), sprintf('%s (%.2f kg)',lbl,m), ...
                        'VerticalAlignment','bottom','Color',clr);
                end
                function addMass(P, field, massVar, tag)
                    if isfield(P,field)
                        plotMass(P.(field), evalin('caller',massVar), [0.25 0.6 0.1], tag);
                    end
                end
        end
        

        function showWingWithChord(obj)
            % plots nodes, element centre-lines AND draws rectangular
            % chord panels so the wing doesn't look like a line
            chord = obj.aeromesh.chord(1,1);
            fem = obj.fem;
            rgParam = obj.rgParam;
            ea = obj.aeromesh.elastic_axis(1,1);

            clf,  hold on, grid on, axis equal
            view(40,20)
            xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
            title('Pazy wing mesh with aerodynamic chord')

            XYZ = obj.ensureXYZ();      % [nNode×3]
            N   = obj.fem.num_node;

            % === (a) draw nodes (blue) =================================
            scatter3(XYZ(:,1),XYZ(:,2),XYZ(:,3),25,'b','filled')

            % === (b) draw elastic-axis beam (black) ====================
            for e = 1:obj.fem.num_elem
                idx = obj.fem.connectivities(e,[1 3 2]); % {0,2,1} order
                plot3(XYZ(idx,1),XYZ(idx,2),XYZ(idx,3),'-k','LineWidth',1.2)
            end

            % === (c) draw chord rectangles =============================
            %  local chord runs +x  (because Pazy beam is laid out in –y
            %  span-wise direction).
            ex = [1 0 0];                                  % chord dir (A-frame)
            ea2le =  ea*chord;              % EA → LE
            ea2te = (1-ea)*chord;           % EA → TE

            for n = 1:N-1        % panel between node-n and node-n+1
                r1 = XYZ(n  ,:);
                r2 = XYZ(n+1,:);
                % Four corner pts of rectangle
                P1 = r1 - ea2le*ex;   % leading-edge, inboard
                P2 = r2 - ea2le*ex;   % leading-edge, outboard
                P3 = r2 + ea2te*ex;   % trailing
                P4 = r1 + ea2te*ex;
                patch('XData',[P1(1) P2(1) P3(1) P4(1)], ...
                      'YData',[P1(2) P2(2) P3(2) P4(2)], ...
                      'ZData',[P1(3) P2(3) P3(3) P4(3)], ...
                      'FaceColor',obj.colWing, 'FaceAlpha',0.35, ...
                      'EdgeColor','none');
            end
                     % remove outline for clarity
            af = obj.aeromesh.airfoils.one;        % 2×N, (x/c,z/c)
    
            for n = 1:fem.num_node
                % Elastic-axis nodes lie on Y≈0 (because Pazy is unswept).
                % if abs(fem.coordinates(n,2)) < 1e-3
                    clrs =  [0.2 0.6 1];
                    obj.addAirfoilPatch(af, fem.coordinates(n,:), chord,clrs, 0.30);
                % end
            end
        end

        %--------------------------------------------------------------
        function showWingAndRigidBody(obj)
            % calls wing plot then overlays lumped masses
            obj.showWingWithChord();
            ax = gca;

            prm = obj.rgParam;
            % Plot the fuselage as a cylinder centered at the wing root
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
            d2 = set(h2, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.9);  % red sphere
            % legend
            % Tail mass
            tailPos = obj.paramsRB.tail.pos;
            h3 = surf(ax, sx + tailPos(1), sy + tailPos(2), sz + tailPos(3));
            d3 = set(h3, 'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.9);  % green sphere
            % legend('fuselage CM', h1, 'battery CM', h2, 'Tail CM', h3 )
            % legend(h1, 'fuselage CM')
            % legend( d2,'Battery CM' )
            % legend( d3,'Tail CM' )
            % legend('h1', 'h2','h3' )
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
            % --- hide all surfaces from the legend (fuselage tube + spheres)
            set(findall(ax,'Type','Surface'),'HandleVisibility','off');        % generic
            % (or explicitly:)
            % fusSurf.HandleVisibility = 'off';  % if you kept the fuselage surf handle
            % set([h1 h2 h3 h4 h5],'HandleVisibility','off');
            
            % --- legend proxy markers (circles), no lines
            hold(ax,'on');
            mSize = 80;   % legend marker size
            pFus  = scatter3(ax, NaN,NaN,NaN, mSize, 'o', 'filled', ...
                             'MarkerFaceColor','k','MarkerEdgeColor','k', ...
                             'DisplayName','Fuselage CM');
            pTail = scatter3(ax, NaN,NaN,NaN, mSize, 'o', 'filled', ...
                             'MarkerFaceColor','g','MarkerEdgeColor','k', ...
                             'DisplayName','Tail CM');
            pBat  = scatter3(ax, NaN,NaN,NaN, mSize, 'o', 'filled', ...
                             'MarkerFaceColor','r','MarkerEdgeColor','k', ...
                             'DisplayName','Battery CM');

            pWing = scatter3(ax, NaN,NaN,NaN, mSize, 'o', 'filled', ...
                             'MarkerFaceColor',[1 0 1],'MarkerEdgeColor','k', ...
                             'DisplayName','Wing Ctrl Surfaces');
            
            % --- build the legend using the proxies (shows circles, not rectangles)
            lgd = legend(ax, [pFus pTail pBat pWing], {'Fuselage CM','Tail CM','Battery CM','Wing Ctrl Surfaces'}, ...
                         'Location','best','Interpreter','none');
            lgd.Box = 'on';
            % Finalize axes
            view(30, 15);  % set a good 3D view angle (azimuth, elevation)
            grid on;
            hold(ax, 'off');
            % % All mass points ------------------------------------------------
            % masses = { ...
            %    'Fuselage', prm.mass*0 + prm.mass, prm.rCM.' ;  % show CM
            %    'Tail'    , 0.12 ,  [-0.70  0  0] ;
            %    'Battery' , 0.18 ,  [-0.10  0  0] ;
            %    'R-tip'   , 0.035, [ 0.60  0.55 0] ;
            %    'L-tip'   , 0.035, [ 0.60 -0.55 0] };
            % 
            % for k = 1:size(masses,1)
            %     m  = masses{k,2};
            %     r  = masses{k,3};
            %     [xs,ys,zs] = sphere(8);
            %     sc = 0.03*(m/0.1)^(1/3);   % scale marker ~ m^(1/3)
            %     surf(sc*xs + r(1), sc*ys + r(2), sc*zs + r(3), ...
            %          'FaceColor',[0.9 0.3 0.1], 'EdgeColor','none', ...
            %          'FaceAlpha',0.7)
            %     text(r(1),r(2),r(3)+0.03,sprintf('%s',masses{k,1}), ...
            %          'FontSize',8,'Color','k');
            % end
            % legend({'beam nodes','elastic axis','wing surf','masses'}, ...
            %        'Location','bestoutside')
            % plot3(prm.geom.fuselageLine(1,:), ...
            % prm.geom.fuselageLine(2,:), ...
            % prm.geom.fuselageLine(3,:), 'k','LineWidth',2);
        end

        %--------------------------------------------------------------
        function XYZ = ensureXYZ(obj)
            XYZ = obj.fem.coordinates;
            % -> row-wise [n×3]
            if size(XYZ,2)~=3, XYZ = XYZ.'; end
        end
        function addAirfoilPatch(obj, afXY, nodePosA, chord, clr, alphaFace)
            % addAirfoilPatch  Draws a semi-transparent airfoil at a given node position.
            %
            % INPUTS
            %   afXY        2×N or N×2   – air-foil coordinates in (x/c ,  z/c ).
            %   nodePosA    1×3          – node coordinates in A-frame [m].
            %   chord       scalar [m]   – local chord length.
            %   clr         1×3 RGB      – patch colour (e.g. [0.3 0.6 1]).
            %   alphaFace   scalar 0-1   – patch face transparency.
            %
            % USAGE EXAMPLE  (inside your plot routine):
            %   for k = 1:fem.num_node
            %       if abs(fem.coordinates(k,2)) < 1e-3        % nodes on elastic axis
            %           addAirfoilPatch(af, fem.coordinates(k,:), 0.10, [0.2 0.6 1], .25);
            %       end
            %   end
            %
            % -------------------------------------------------------------------------
            
            % --- ensure afXY is 2×N -----------------------------------------------
            if size(afXY,1) ~= 2
                afXY = afXY.';                           % transpose if needed
            end
            x_local = afXY(1,:)*chord;                   % scale to chord
            z_local = afXY(2,:)*chord;
            
            % airfoil lives in the X–Z plane of body-A frame (positive X downstream).
            X = nodePosA(1) + x_local;
            Y = nodePosA(2)*ones(size(X));               % keep the spanwise coord
            Z = nodePosA(3) + z_local;
            
            patch(X, Y, Z, 1, ...
                  'FaceColor', clr, ...
                  'FaceAlpha', alphaFace, ...
                  'EdgeColor', 'none');                  % remove outline for clarity
        end

        function activeDofs = defineActiveDofs_pazyWing(fem, dofBC)
        %DEFINEACTIVEDOFS_PAZYWING  Identify the master (active) DOFs for condensation
        %
        % Inputs:
        %   numNodes    = total number of FE nodes in the Pazy wing mesh
        %   dofsPerNode = number of DOFs per node (6 if 3D beam, 3 if 2D)
        %   clampedNode = ID of the clamped (root) node (often node 0)
        %   lumpNodes   = array of node IDs that have lumped masses or rods
        %                 e.g. [5, 9, 12] if your model lumps are placed there
        %
        % Output:
        %   activeDofs  = column vector listing the global DOFs that remain
        %                 after condensation. (All others are 'omitted' or slaved.)
            % disp(fem)
            activeDofs = [];
            disp(fem)
            dofPerNode = 6;
            
            for n = 1 : fem.num_node
                % If node is clamped, skip it entirely
                if ismember(n, dofBC)
                    continue
                end
                
                % If node is among the lumps, keep all its DOFs
                if ismember(n, fem.lumped_mass_nodes)
                    baseIndex = double((n-1)*dofPerNode); % (1-based "n", so subtract 1)
                    theseDofs = baseIndex + (1 : dofPerNode);
                    activeDofs = [activeDofs, theseDofs]; %#ok<AGROW>
                else
                    % % Possibly keep only translations or skip entirely,
                    % % % e.g. keep only translation DOFs for non-lump nodes:
                    %   baseIndex = double((n-1)*dofPerNode);
                    %   transDofs = baseIndex + [1 2 3]; % e.g. x,y,z only
                    %   activeDofs = [activeDofs, transDofs];
                    % 
                    % Or skip them altogether:
                    % do nothing
                end
            end
        
            % remove duplicates & sort:
            activeDofs = unique(activeDofs, 'sorted');
        
        end

    end
end

