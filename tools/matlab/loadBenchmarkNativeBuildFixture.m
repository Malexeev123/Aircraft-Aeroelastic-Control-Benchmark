function [fixture,fixturePath] = loadBenchmarkNativeBuildFixture(kind,repositoryRoot)
%LOADBENCHMARKNATIVEBUILDFIXTURE Load a hash-locked MEX build fixture.
%   The compact fixtures preserve the accepted code-generation dimensions and
%   parity inputs without distributing ordinary simulation checkpoints.

arguments
    kind (1,1) string {mustBeMember(kind,["fixed_interval","scheduled"])}
    repositoryRoot (1,1) string
end

switch kind
    case "fixed_interval"
        name = "phase18c_v17a_fixed_interval_build_fixture_v1.mat";
        expectedFileHash = ...
            "7e7dcb82bca0aa3d867f880b100c6e6b45e435c2b4ad95f6271eab9c6858a648";
        expectedSchema = "pazy-native-fixed-interval-build-fixture-v1";
        expectedCheckpointHash = ...
            "21b8970da975a0d7c469923d57a2f293abedf17b932a33700c0ff1c30f981f54";
        required = ["provider","state","disturbance","context","packet"];
        expectedStateCount = 121;
    case "scheduled"
        name = "phase18c_v17a_scheduled_build_fixture_v1.mat";
        expectedFileHash = ...
            "4d32586a5ccee22cc5c166e1c332a6bd19b4e1b83e91fa22328e15925779d2b5";
        expectedSchema = "pazy-native-scheduled-build-fixture-v1";
        expectedCheckpointHash = ...
            "daafa0cd9b144180fcab2ce697e716ee6bd0844a024885a19505fda5bba2107b";
        required = ["model","packet","context","state"];
        expectedStateCount = 620;
end

fixturePath = fullfile(repositoryRoot,"MatlabFlex","configs", ...
    "benchmark","native-build-fixtures",name);
assert(isfile(fixturePath),"AeroFlex:NativeBuildFixtureMissing", ...
    "The required native-build fixture is missing: %s",fixturePath);
actualHash = localFileHash(fixturePath);
assert(actualHash==expectedFileHash,"AeroFlex:NativeBuildFixtureHash", ...
    "The native-build fixture hash does not match the release contract.");

% Copy before loading because Windows MATLAB can intermittently reject a
% direct MAT-file read through a WSL UNC path. The source bytes remain fixed.
temporary = string(tempname)+".mat";
[copied,message] = copyfile(fixturePath,temporary,"f");
assert(copied,"AeroFlex:NativeBuildFixtureCopy", ...
    "Unable to copy the native-build fixture: %s",message);
cleanup = onCleanup(@()localDelete(temporary));
loaded = load(temporary,"fixture");
assert(isfield(loaded,"fixture") && isstruct(loaded.fixture), ...
    "AeroFlex:NativeBuildFixturePayload", ...
    "The native-build fixture payload is missing or invalid.");
fixture = loaded.fixture;
clear cleanup
localDelete(temporary)

assert(isfield(fixture,"schemaVersion") && ...
    string(fixture.schemaVersion)==expectedSchema, ...
    "AeroFlex:NativeBuildFixtureSchema", ...
    "The native-build fixture schema does not match its declared role.");
assert(isfield(fixture,"sourceCheckpointSha256") && ...
    string(fixture.sourceCheckpointSha256)==expectedCheckpointHash, ...
    "AeroFlex:NativeBuildFixtureProvenance", ...
    "The native-build fixture source-checkpoint hash is not accepted.");
assert(all(isfield(fixture,required)), ...
    "AeroFlex:NativeBuildFixtureFields", ...
    "The native-build fixture does not contain every required field.");
assert(isnumeric(fixture.state) && iscolumn(fixture.state) && ...
    numel(fixture.state)==expectedStateCount && ...
    all(isfinite(fixture.state)), ...
    "AeroFlex:NativeBuildFixtureState", ...
    "The native-build fixture state is dimensionally invalid or nonfinite.");
end

function value = localFileHash(path)
file = fopen(path,"rb");
assert(file>=0,"AeroFlex:NativeBuildFixtureRead", ...
    "Cannot read the native-build fixture: %s",path);
cleanup = onCleanup(@()fclose(file));
engine = java.security.MessageDigest.getInstance("SHA-256");
while ~feof(file)
    bytes = fread(file,1024*1024,"*uint8");
    if isempty(bytes), break, end
    engine.update(bytes);
end
value = lower(string(sprintf("%02x",typecast( ...
    engine.digest(),"uint8"))));
clear cleanup
end

function localDelete(path)
if isfile(path), delete(path); end
end
