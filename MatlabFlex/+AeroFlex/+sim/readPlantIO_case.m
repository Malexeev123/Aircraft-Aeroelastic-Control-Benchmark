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
            % Temporary fallback if your current PlantRunTime only has
            % readSensors:
            [z_k, t_k] = plant.readSensors(cfg, Ssim.y);

            io = struct();
            io.t     = t_k;
            io.yWing = z_k;
            io.yRB   = [];      % TODO: populate from rigid-body sensors
            io.ref   = [];      % TODO: guidance reference

        otherwise
            error('Unknown bodyCase = "%s".', bodyCase);
    end
end