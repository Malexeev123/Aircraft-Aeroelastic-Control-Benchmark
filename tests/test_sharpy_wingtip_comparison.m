function tests = test_sharpy_wingtip_comparison
%TEST_SHARPY_WINGTIP_COMPARISON Verify the supplied comparison dataset.
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
project = setupProject(ValidateEntryPoints=true);
root = string(project.repositoryRoot);
testCase.TestData.dataPath = fullfile(root,"results","validation", ...
    "sharpy-wingtip-comparison","data","wingtip_comparison.csv");
testCase.TestData.summaryPath = fullfile(root,"results","validation", ...
    "sharpy-wingtip-comparison","wingtip_comparison_summary.json");
end


function testAcceptedComparisonMetrics(testCase)
data = readtable(testCase.TestData.dataPath, ...
    "VariableNamingRule","preserve");
summary = jsondecode(fileread(testCase.TestData.summaryPath));

sharpy = double(data.sharpy_mean_tip_m(:));
matlab = double(data.matlab_symmetric_tip_m(:));
error = matlab-sharpy;

verifyEqual(testCase,height(data),summary.sampleCount);
verifyTrue(testCase,all(diff(data.time_s)>0));
verifyTrue(testCase,all(isfinite(data{:,vartype("numeric")}),"all"));
verifyEqual(testCase,max(sharpy),summary.metrics.sharpyPeakMeters, ...
    "AbsTol",summary.integrityTolerance.metricAbsolute);
verifyEqual(testCase,max(matlab),summary.metrics.matlabPeakMeters, ...
    "AbsTol",summary.integrityTolerance.metricAbsolute);
verifyEqual(testCase,max(matlab)/max(sharpy), ...
    summary.metrics.peakRatioMatlabToSharpy, ...
    "AbsTol",summary.integrityTolerance.ratioAbsolute);
verifyEqual(testCase,sqrt(mean(error.^2)), ...
    summary.metrics.rmsErrorMeters, ...
    "AbsTol",summary.integrityTolerance.metricAbsolute);
verifyEqual(testCase,max(abs(error)), ...
    summary.metrics.peakAbsoluteErrorMeters, ...
    "AbsTol",summary.integrityTolerance.metricAbsolute);
end


function testComparisonContractIsTrimRelative(testCase)
summary = jsondecode(fileread(testCase.TestData.summaryPath));

verifyEqual(testCase,string(summary.comparisonFrame), ...
    "initial_trim_relative_tip_displacement");
verifyEqual(testCase,string(summary.projectionPolicy), ...
    "disabled_for_wing_only_replay");
verifyEqual(testCase,summary.flightCondition.speedMetersPerSecond,40);
verifyEqual(testCase,summary.flightCondition.angleOfAttackDegrees,1);
verifyEqual(testCase,summary.gust.peakMetersPerSecond,0.8);
end
