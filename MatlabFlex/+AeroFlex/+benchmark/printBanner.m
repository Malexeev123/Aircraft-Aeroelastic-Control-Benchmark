function printBanner(plan,stage,details)
%PRINTBANNER Print a compact benchmark startup or completion summary.

arguments
    plan (1,1) struct
    stage (1,1) string {mustBeMember(stage,["start","complete"])}
    details (1,1) struct = struct()
end

line = repmat('=',1,70);
fprintf("\n%s\n",line);
if stage=="start"
    fprintf(" PAZY AEROELASTIC CONTROL BENCHMARK\n");
    fprintf(" Coupled nonlinear aeroelastic estimation and control\n");
    fprintf(" Model: %-6s  Case: %-4s  Status: %s\n", ...
        plan.modelRevision,plan.caseId,plan.qualificationStatus);
    fprintf(" Controller: %-10s  Gust: %-8s  Duration: %.3f s\n", ...
        plan.controllerMode,localOnOff(plan.gustEnabled), ...
        plan.durationSeconds);
    fprintf(" Native kernels: %-8s  Diagnostics: level %d\n", ...
        plan.nativeKernelPolicy,plan.diagnosticsLevel);
    fprintf(" Output root: %s\n",plan.outputRoot);
    if ~plan.qualified
        fprintf(" WARNING: This case is not a closed benchmark member.\n");
        fprintf("          %s\n",plan.unavailableReason);
    end
else
    fprintf(" PAZY BENCHMARK RUN COMPLETE\n");
    fprintf(" Case: %-4s  Qualification: %s\n", ...
        plan.caseId,plan.qualificationStatus);
    if isfield(details,"status")
        fprintf(" Result status: %s\n",string(details.status));
    end
    if isfield(details,"onlineWallSeconds")
        fprintf(" Online simulation time: %.3f s\n", ...
            double(details.onlineWallSeconds));
    end
    if isfield(details,"runnerExecutionWallSeconds")
        fprintf(" Runner execution time: %.3f s\n", ...
            double(details.runnerExecutionWallSeconds));
    end
    if isfield(details,"preparationSeconds")
        fprintf(" Preparation time: %.3f s (excluded from online timing)\n", ...
            double(details.preparationSeconds));
    end
    if isfield(details,"totalEntryWallSeconds")
        fprintf(" Total entry time: %.3f s\n", ...
            double(details.totalEntryWallSeconds));
    end
    if isfield(details,"runRoot")
        fprintf(" Results: %s\n",string(details.runRoot));
    end
end
fprintf("%s\n",line);
end

function value = localOnOff(flag)
if flag
    value = "enabled";
else
    value = "disabled";
end
end
