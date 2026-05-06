classdef TestTrim < matlab.unittest.TestCase
%DETERMINISTIC TRIM TEST  (no wind, steady flight)

    methods (Test)
        function steadyTrim(testCase)
            cfg = nominalConfig;
            sim = SimRunner(cfg);
            sim.run();
            x = sim.log.x;
            dx = max(abs(diff(x,1,2)),[],'all');
            testCase.verifyLessThan(dx,1e-10,...
                'State should be constant in deterministic trim.');
        end
    end
end

