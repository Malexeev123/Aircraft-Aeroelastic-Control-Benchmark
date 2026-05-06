classdef TestnMPCStep < matlab.unittest.TestCase
%NMPC STEP–RESPONSE TO COMMANDED PITCH

    methods (Test)
        function stepResponse(testCase)
            cfg = nominalConfig;
            cfg.controller = nMPC(cfg.mpc);
            sim = SimRunner(cfg);
            % change reference at half‑time
            tmid = cfg.t_end/2;
            cfg.mpc.refFcn = @(t) (t>tmid)*deg2rad(5);   % step 5°
            sim.run();

            pitch = sim.log.x(cfg.idxPitch,:);
            overshoot = max(pitch) - deg2rad(5);
            testCase.verifyLessThan(overshoot,deg2rad(1),...
                'nMPC overshoot should be <1°.');
        end
    end
end

