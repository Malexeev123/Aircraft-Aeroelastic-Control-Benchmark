function sm_logging(onOff)
logsout = Simulink.SimulationData.Dataset; %#ok<NASGU>
if strcmpi(onOff,'on'), set_param(bdroot,'SignalLogging','on');
else,                   set_param(bdroot,'SignalLogging','off'); end
end
