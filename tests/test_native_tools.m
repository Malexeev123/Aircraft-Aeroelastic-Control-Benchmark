function tests = test_native_tools
%TEST_NATIVE_TOOLS Build-contract and numerical-parity tests for native tools.
%   Run with:
%       results = runtests("tests/test_native_tools.m");
%       assertSuccess(results)
%
%   Missing kernels are built from the public, hash-locked fixtures. Existing
%   compatible caches are reused. Every builder compares its C++ MEX output
%   with the corresponding MATLAB implementation before reporting success.

tests = functiontests(localfunctions);
end


function setupOnce(testCase)
repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repositoryRoot);

project = setupProject(ValidateEntryPoints=true, ...
    ChangeCurrentFolder=false);

report = buildBenchmarkTools(Force=false,RunParity=true, ...
    TimingRepetitions=1,ProjectInfo=project);

testCase.TestData.project = project;
testCase.TestData.report = report;
end


function testSupportedNativeEnvironment(testCase)
environment = testCase.TestData.report.environment;

verifyTrue(testCase,environment.releaseSupported, ...
    "Native kernels are qualified for MATLAB R2025b.");
verifyTrue(testCase,environment.architectureSupported, ...
    "Native kernels are qualified for 64-bit Windows MATLAB.");
verifyTrue(testCase,environment.matlabCoderAvailable, ...
    "MATLAB Coder is required to build the native kernels.");
verifyTrue(testCase,environment.cppCompilerSelected, ...
    "Select a supported compiler with 'mex -setup C++'.");
end


function testAllNativeFamiliesPassParity(testCase)
report = testCase.TestData.report;
expected = [ ...
    "fixed_interval"
    "scheduled_interval"
    "scheduled_horizon"
    "scheduled_value_horizon"
    "scheduled_causal_rollout"];

verifyTrue(testCase,report.passed);
verifyEqual(testCase,string({report.components.name}).',expected);
verifyTrue(testCase,all([report.components.passed]));

for index = 1:numel(report.components)
    details = report.components(index).details;
    verifyTrue(testCase,isfield(details,"passed"));
    verifyTrue(testCase,logical(details.passed), ...
        sprintf("Native family did not pass: %s",expected(index)));
end
end


function testCacheManifestsRemainActive(testCase)
project = setupProject(ValidateEntryPoints=true, ...
    ChangeCurrentFolder=false);
check = buildBenchmarkTools(Action="check",ProjectInfo=project);

verifyTrue(testCase,check.passed);
verifyEqual(testCase,string({check.components.name}).', ...
    string({testCase.TestData.report.components.name}).');
verifyTrue(testCase,all([check.components.active]));
verifyTrue(testCase,all(strlength(string( ...
    {check.components.path})) > 0));
end


function testNativeRequiredInstallationPasses(testCase)
project = setupProject(ValidateEntryPoints=true, ...
    ChangeCurrentFolder=false);
status = verifyBenchmarkInstallation(RequireNativeKernels=true, ...
    PrintSummary=false,ProjectInfo=project);

verifyTrue(testCase,status.passed);
verifyTrue(testCase,status.nativeRequired);
verifyTrue(testCase,status.nativeTools.passed);
verifyTrue(testCase,status.registry.passed);
end
