classdef TestEKF < matlab.unittest.TestCase
%KF VS TRUE STATES UNDER SENSOR NOISE

    methods (Test)
        function ekfConverges(testCase)
            cfg = nominalConfig;
            cfg.estimator = EKF(cfg.kf);               % inject KF
            sim = SimRunner(cfg);
            sim.run();
            err = sim.log.x - sim.estimator.xhat;
            rms = rms(err,'all');
            testCase.verifyLessThan(rms,0.05,...
                'EKF should converge to <5 % RMS error.');
        end
    end
end
