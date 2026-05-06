
classdef PID
    properties
        Kp Ki Kd
    end
    methods
        function obj = PID(kp,ki,kd)
            obj.Kp=kp; obj.Ki=ki; obj.Kd=kd;
        end
        function u = computeInput(obj,~,error,derror)
            u = obj.Kp*error + obj.Ki*0 + obj.Kd*derror;
        end
    end
end
