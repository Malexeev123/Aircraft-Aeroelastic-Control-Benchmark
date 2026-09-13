function tests = test_benchmark_runner_dispatch
%TEST_BENCHMARK_RUNNER_DISPATCH Check public plans reach their runtime adapters.
tests = functiontests(localfunctions);
end

function setupOnce(~)
setupProject(PrintSummary=false);
end

function testReadmeCustomScenarioReachesAdapter(testCase)
reference = struct("timeSeconds",[0,5,10],"speedMps",[15,17,20], ...
    "pitchRad",deg2rad([10,10.5,10]), ...
    "phaseLabel",["attitude","combined","speed"]);
definition = AeroFlex.benchmark.customCaseDefinition( ...
    Name="combined_maneuver",DurationSeconds=10,Reference=reference);
[summary,plan] = runCustomBenchmarkCase(definition,Execute=false);
verifyEqual(testCase,summary.status,"PLAN_ONLY");
verifyEqual(testCase,plan.runnerKind,"custom_scheduled");
args = AeroFlex.benchmark.internal.runnerArguments(plan,"unused_test_output");
verifyEqual(testCase,args{1},10);
verifyEqual(testCase,args{end},definition);
verifyEqual(testCase,args{43},definition.gust.enabled);
% Stop at the real runner's argument validation, before any simulation.
plan.durationSeconds = -1;
verifyError(testCase,@() AeroFlex.benchmark.internal.executePlan( ...
    plan,"unused_test_output"),'MATLAB:validators:mustBePositive');
end

function testFormalScheduledMembersReachAdapter(testCase)
for caseId = ["B1","B2"]
    [~,plan] = runBenchmarkCase(caseId,Execute=false);
    args = AeroFlex.benchmark.internal.runnerArguments(plan,"unused_test_output");
    verifyEqual(testCase,plan.runnerKind,"scheduled_case_b");
    verifyEqual(testCase,args{1},plan.durationSeconds);
    verifyEqual(testCase,args{43},plan.gustEnabled);
    verifyEmpty(testCase,fieldnames(args{end}));
    plan.durationSeconds = -1;
    verifyError(testCase,@() AeroFlex.benchmark.internal.executePlan( ...
        plan,"unused_test_output"),'MATLAB:validators:mustBePositive');
end
end

function testCaseAArgumentsRemainUnchanged(testCase)
for caseId = ["A1","A2","A3"]
    [~,plan] = runBenchmarkCase(caseId,Execute=false);
    args = AeroFlex.benchmark.internal.runnerArguments(plan,"unused_test_output");
    verifyEqual(testCase,plan.runnerKind,"formal_case_a");
    verifyEqual(testCase,numel(args),12);
    verifyEqual(testCase,args{1},plan.durationSeconds);
    verifyEqual(testCase,args{2},plan.subcase);
    verifyEqual(testCase,args{end},plan.runtimeAcceleration);
end
end

function testDocumentedModelFamiliesPlan(testCase)
for bodyCase = ["wingOnly","coupledFull"]
    for mode = ["openloop","nmhe_nmpc"]
        [result,plan] = runPazyModelWorkflow(BodyCase=bodyCase, ...
            SimulationMode=mode,DurationSeconds=0.1,Execute=false);
        verifyEqual(testCase,result.status,"PLAN_ONLY");
        verifyEqual(testCase,plan.bodyCase,bodyCase);
        verifyEqual(testCase,plan.simulationMode,mode);
    end
end
end

function testUnsupportedCustomScopeFailsBeforeDispatch(testCase)
definition = AeroFlex.benchmark.customCaseDefinition( ...
    Reference=struct("timeSeconds",[0,5],"speedMps",[15,20],"altitudeM",[0,1]));
[~,plan] = runCustomBenchmarkCase(definition,Execute=false);
verifyFalse(testCase,plan.executionAllowed);
verifyError(testCase,@() runCustomBenchmarkCase(definition,Execute=true), ...
    'AeroFlex:CustomCaseExecutionUnavailable');
end

function testUnknownRunnerStillRejected(testCase)
plan = struct('runnerKind',"unknown",'caseId',"unknown");
verifyError(testCase,@() AeroFlex.benchmark.internal.runnerArguments( ...
    plan,"unused_test_output"),'AeroFlex:BenchmarkCaseUnavailable');
end
