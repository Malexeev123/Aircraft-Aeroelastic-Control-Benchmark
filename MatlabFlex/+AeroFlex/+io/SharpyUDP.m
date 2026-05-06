
classdef SharpyUDP
    properties
        host
        port
        udpObj
    end
    methods
        function obj = SharpyUDP(host,port)
            if nargin<1, host='127.0.0.1'; end
            if nargin<2, port=9000; end
            obj.host = host; obj.port = port;
            % obj.udpObj = udpport('byte', 'IPV4'); % R2020b+
        end
        function send(obj,data)
            %#ok
        end
        function data = recv(obj)
            data = [];
        end
    end
end
