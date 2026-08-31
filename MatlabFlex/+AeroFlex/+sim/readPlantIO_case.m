function io = readPlantIO_case(plant, cfg, Ssim, bodyCase)
%READPLANTIO_CASE Unified sensor interface for wing-only and coupled cases.

    switch lower(string(bodyCase))

        case "wingonly"
            [z_k, t_k] = plant.readSensors(cfg, Ssim.y);

            io = struct();
            io.t     = t_k;
            io.yWing = z_k;
            io.yRB   = [];
            io.ref   = [];

        case "coupledfull"
            %   API:
            %
            %   io = plant.readCoupledSensors(cfg, Ssim.y);
            %
            % where io contains:
            %   io.t
            %   io.yWing
            %   io.yRB
            %   io.rbState
            %   io.alpha
            %   io.Uinf
            %
            % Coupled propagation advances xFlex, whereas x is retained for
            % the wing-only path. A real measurement must use the propagated
            % coupled q1 coordinates; retained logged sensing uses its
            % existing PlantRunTime fallback.
            if cfg.forceRealSense
                q1 = plant.xFlex(plant.model.idx.q1);
                z_k = plant.sensor.measure(q1);
                t_k = plant.t;
            else
                [z_k, t_k] = plant.readSensors(cfg, Ssim.y);
            end

            io = struct();
            io.t     = t_k;
            io.yWing = z_k;
            io.yRB   = [];      % TODO: populate from rigid-body sensors
            io.ref   = [];      % TODO: guidance reference

        otherwise
            error('Unknown bodyCase = "%s".', bodyCase);
    end
end
