function tests = test_general_model_workflows
%TEST_GENERAL_MODEL_WORKFLOWS Clean-clone model-data and entry-point checks.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
project = setupProject(ValidateEntryPoints=true);
testCase.TestData.root = string(project.repositoryRoot);
testCase.TestData.project = project;
end

function testCuratedModelAssetsMatchManifest(testCase)
status = AeroFlex.benchmark.verifyGeneralModelAssets( ...
    testCase.TestData.root);
verifyTrue(testCase,status.passed);
verifyEqual(testCase,status.assetCount,8);
verifyTrue(testCase,all([status.records.exists]));
verifyTrue(testCase,all([status.records.passed]));
end

function testMissingAssetIsDetected(testCase)
temporaryRoot = string(tempname);
mkdir(fullfile(temporaryRoot,"MatlabFlex","configs","benchmark"));
cleanup = onCleanup(@()rmdir(temporaryRoot,"s"));
sourceManifest = fullfile(testCase.TestData.root,"MatlabFlex", ...
    "configs","benchmark","pazy_general_model_assets_v1.json");
targetManifest = fullfile(temporaryRoot,"MatlabFlex", ...
    "configs","benchmark","pazy_general_model_assets_v1.json");
copyfile(sourceManifest,targetManifest);

status = AeroFlex.benchmark.verifyGeneralModelAssets(temporaryRoot);
verifyFalse(testCase,status.passed);
verifyFalse(testCase,any([status.records.exists]));
clear cleanup
end

function testDynamicRuntimeAssetsIncludePhysicalChartContracts(testCase)
status = prepareBenchmarkReleaseAssets(Action="check", ...
    ProjectInfo=testCase.TestData.project,PrintSummary=false);
roles = string({status.records.role}).';
selected = startsWith(roles,"v17a_source_selectedContract:");

verifyTrue(testCase,status.passed);
verifyEqual(testCase,sum(selected),29);
verifyTrue(testCase,all([status.records(selected).exists]));
verifyTrue(testCase,all([status.records(selected).passed]));
contractPaths = string({status.records(selected).path}).';
verifyTrue(testCase,all(contains(contractPaths, ...
    "fixed_node27_contract") & endsWith(contractPaths,".h5")));
end

function testSharedWorkflowPlans(testCase)
requests = { ...
    {"wingOnly","openloop",true}; ...
    {"wingOnly","nmhe_nmpc",true}; ...
    {"coupledFull","openloop",false}};
for index = 1:numel(requests)
    request = requests{index};
    [result,plan] = runPazyModelWorkflow( ...
        BodyCase=request{1},SimulationMode=request{2}, ...
        GustEnabled=request{3},DurationSeconds=0.05,Execute=false);
    verifyEqual(testCase,result.status,"PLAN_ONLY");
    verifyFalse(testCase,plan.executed);
    verifyEqual(testCase,plan.durationSeconds,0.05,"AbsTol",0);
    if request{1}=="coupledFull"
        verifyEqual(testCase,plan.operatingPointPolicy, ...
            "production_default");
        verifyEqual(testCase,plan.libraryMode,"production");
        verifyEqual(testCase,plan.initializationMode, ...
            "qualified_source_trim");
        verifyEqual(testCase,plan.targetSpeedMps,15,"AbsTol",0);
        verifyEqual(testCase,plan.targetAngleOfAttackDeg,10,"AbsTol",0);
    else
        verifyEqual(testCase,plan.operatingPointPolicy,"supplied_source");
        verifyEqual(testCase,plan.libraryMode,"supplied_source");
        verifyEqual(testCase,plan.initializationMode,"recompute_trim");
    end
end
end

function testExplicitOperatingPointPlan(testCase)
[~,plan] = runPazyModelWorkflow(BodyCase="coupledFull", ...
    SimulationMode="openloop",GustEnabled=false,Execute=false, ...
    TargetSpeedMps=40,TargetAngleOfAttackDeg=1);
verifyEqual(testCase,plan.operatingPointPolicy,"explicit");
verifyEqual(testCase,plan.targetSpeedMps,40,"AbsTol",0);
verifyEqual(testCase,plan.targetAngleOfAttackDeg,1,"AbsTol",0);
verifyEqual(testCase,plan.libraryMode,"production");
verifyEqual(testCase,plan.initializationMode,"qualified_source_trim");
end

function testProductionInitialStateUsesQualifiedRuntimeChart(testCase)
[library,source] = AeroFlex.benchmark.loadProductionScheduledLibrary( ...
    testCase.TestData.root,[15,10]);
expected = AeroFlex.sched.evalLibrary(library,[15,10],struct( ...
    "requireCompatible",true,"noExtrapolate",true, ...
    "interpTol",1e-10,"fullCoordinateRuntimeCandidate", ...
    source.fullCoordinateRuntimeCandidate));
sourceIndex = find(max(abs(library.mu-[15,10]),[],2)<=1e-10);

verifyTrue(testCase,isscalar(sourceIndex));
verifyEqual(testCase,source.p5.pointIds,sourceIndex);
verifyEqual(testCase,source.p5.x_eq,expected.x_eq,"AbsTol",0);
verifyEqual(testCase,source.trim.states,expected.x_eq(:),"AbsTol",0);
verifyGreaterThan(testCase,norm( ...
    library.points(sourceIndex).x_eq(:)-expected.x_eq(:),inf),0.5);
verifyGreaterThanOrEqual(testCase,source.trim.thrust,0);

for index = 1:numel(library.points)
    artifactPath = string(library.points(index).p5.artifactPath);
    verifyTrue(testCase,isfile(artifactPath));
    verifyTrue(testCase,startsWith(lower(artifactPath), ...
        lower(testCase.TestData.root)));
end

[caseALibrary,caseASource] = ...
    AeroFlex.benchmark.loadProductionScheduledLibrary( ...
    testCase.TestData.root,[40,1]);
caseAIndex = find(max(abs(caseALibrary.mu-[40,1]),[],2)<=1e-10);
verifyTrue(testCase,isscalar(caseAIndex));
verifyFalse(testCase,caseASource.dynamicSchedulingSupported);
verifyEqual(testCase,caseASource.p5.x_eq, ...
    caseALibrary.points(caseAIndex).x_eq,"AbsTol",0);
verifyGreaterThanOrEqual(testCase,caseASource.trim.thrust,0);
verifyTrue(testCase,startsWith(lower(string( ...
    caseASource.p5.p5.artifactPath)),lower(testCase.TestData.root)));
end

function testProductionSetupPreservesClosedLoopTimeGrid(testCase)
runId = "test_closed_loop_grid_" + string(java.util.UUID.randomUUID);
setup = AeroFlex.benchmark.prepareProductionModelSetup( ...
    RepositoryRoot=testCase.TestData.root, ...
    SharpyRoot=fullfile(testCase.TestData.root,"TestBenchPazy"), ...
    BodyCase="wingOnly",SimulationMode="nmhe_nmpc", ...
    GustEnabled=false,DurationSeconds=0.01, ...
    OperatingPoint=[40,1],RunId=runId);
cleanup = onCleanup(@()rmdir(setup.paths.run_dir,"s"));

bundle = load(fullfile(setup.paths.for_matlab,"sim_bundle.mat"),"t","y");
verifySize(testCase,bundle.y,[numel(bundle.t),0]);
verifyEqual(testCase,numel(bundle.t),17);
verifyEqual(testCase,diff(bundle.t), ...
    repmat(setup.cfg.library.frozenPackage.parConst.dt,16,1), ...
    "AbsTol",32*eps(setup.cfg.library.frozenPackage.parConst.dt));
verifyFalse(testCase,setup.cfg.trim.useRateProjection);
verifyTrue(testCase,setup.cfg.library.enable);
verifyEqual(testCase,string(setup.cfg.library.updateMode),"frozenTrim");
verifyFalse(testCase, ...
    setup.cfg.library.fullCoordinateRuntimeCandidate.enabled);
verifyGreaterThanOrEqual(testCase,setup.trim.thrust,0);
verifyTrue(testCase,setup.cfg.sharedModelWorkflow.enabled);
verifyTrue(testCase, ...
    setup.cfg.sharedModelWorkflow.exactPackagePrediction);
verifyEqual(testCase,string( ...
    setup.cfg.sharedModelWorkflow.actuatorCommandFrame), ...
    "trim_relative_increment");
verifyEqual(testCase,string(setup.cfg.library.frozenPackage.name), ...
    string(setup.trim.sourceId));
verifyEqual(testCase,setup.cfg.library.frozenPackage.x_eq(:), ...
    setup.trim.states(:),"AbsTol",0);
clear cleanup
end
