function summary = executePlan(plan,runRoot)
%EXECUTEPLAN Dispatch a resolved plan through its retained numerical owner.

arguments
    plan (1,1) struct
    runRoot (1,1) string
end

switch plan.runnerKind
    case "formal_case_a"
        arguments = AeroFlex.benchmark.internal.runnerArguments(plan,runRoot);
        summary = AeroFlex.benchmark.runtime.casea. ...
            run_phase18c_v17a_reciprocal_formal_casea_member_v1( ...
                arguments{:});

    case {"scheduled_case_b","custom_scheduled"}
        arguments = AeroFlex.benchmark.internal.runnerArguments(plan,runRoot);
        summary = AeroFlex.benchmark.runtime. ...
            run_phase18c_v17a_caseb_integrated_profile_v1(arguments{:});

    otherwise
        error("AeroFlex:BenchmarkCaseUnavailable", ...
            "No execution adapter is available for case %s.",plan.caseId);
end
end
