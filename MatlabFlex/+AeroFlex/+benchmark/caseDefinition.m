function definition = caseDefinition(caseId)
%CASEDEFINITION Return the declared intent and qualification of a case.
%   DEFINITION = AEROFLEX.BENCHMARK.CASEDEFINITION(CASEID) returns the
%   immutable research intent used by the public benchmark interface.

arguments
    caseId (1,1) string
end

requestedCaseId = upper(caseId);
mustBeMember(requestedCaseId, ...
    ["A1","A2","A3","B","B1","B2","C","WING_ONLY"]);
canonicalCaseId = requestedCaseId;
if requestedCaseId=="B"
    canonicalCaseId = "B2";
end
definition = struct( ...
    "schemaVersion","pazy-benchmark-case-definition-v1", ...
    "requestedCaseId",requestedCaseId, ...
    "caseId",canonicalCaseId, ...
    "modelRevision","V17A", ...
    "bodyCase","coupled_full", ...
    "defaultDurationSeconds",5, ...
    "subcase","", ...
    "referenceType","attitude", ...
    "commandEnabled",false, ...
    "gustEnabled",false, ...
    "qualificationStatus","", ...
    "qualified",false, ...
    "executionImplemented",true, ...
    "allowResearchOverride",false, ...
    "limitation","", ...
    "evidence",string.empty(1,0));

caseAEvidence = ...
    "context/audits/phase18c-control-validation/speed/" + ...
    "v17a-condensed-rti-compiled-runtime-owner-promotion-v1/" + ...
    "PHASE18C_V17A_FORMAL_CASEA_CLOSURE_V1.json";
caseBEvidence = ...
    "context/audits/phase18c-control-validation/cases/attempts/" + ...
    "v17a-caseb-state-transport-gate-backtracking-extension-v1/" + ...
    "PHASE18C_V17A_CASEB_V90_QGAM_OWNER_DISPOSITION_V1.json";

switch canonicalCaseId
    case "A1"
        definition.subcase = "command_no_gust";
        definition.commandEnabled = true;
        definition.qualificationStatus = "QUALIFIED_CLOSED";
        definition.qualified = true;
        definition.evidence = caseAEvidence;
    case "A2"
        definition.subcase = "hold_gust";
        definition.gustEnabled = true;
        definition.qualificationStatus = "QUALIFIED_CLOSED";
        definition.qualified = true;
        definition.evidence = caseAEvidence;
    case "A3"
        definition.subcase = "command_gust";
        definition.commandEnabled = true;
        definition.gustEnabled = true;
        definition.qualificationStatus = "QUALIFIED_CLOSED";
        definition.qualified = true;
        definition.evidence = caseAEvidence;
    case "B1"
        definition.defaultDurationSeconds = 31;
        definition.subcase = "speed_change_no_gust";
        definition.referenceType = "speed_altitude_attitude";
        definition.commandEnabled = true;
        definition.qualificationStatus = "EXPERIMENTAL_SCHEDULED_CASE";
        definition.allowResearchOverride = true;
        definition.limitation = ( ...
            "Formal qualification is pending; all physical and numerical " + ...
            "acceptance thresholds remain enforced.");
        definition.evidence = caseBEvidence;
    case "B2"
        definition.defaultDurationSeconds = 31;
        definition.subcase = "speed_change_gust";
        definition.referenceType = "speed_altitude_attitude";
        definition.commandEnabled = true;
        definition.gustEnabled = true;
        definition.qualificationStatus = "EXPERIMENTAL_SCHEDULED_CASE";
        definition.allowResearchOverride = true;
        definition.limitation = ( ...
            "Formal qualification is pending; the frozen gust and all " + ...
            "physical and numerical acceptance thresholds remain enforced.");
        definition.evidence = caseBEvidence;
    case "C"
        definition.defaultDurationSeconds = 31;
        definition.subcase = "longitudinal_trajectory_gust";
        definition.referenceType = "longitudinal_trajectory";
        definition.commandEnabled = true;
        definition.gustEnabled = true;
        definition.executionImplemented = false;
        definition.qualificationStatus = "UNAVAILABLE_NOT_FROZEN";
        definition.limitation = ...
            "Formal Case C is deferred until its maneuver and gates are frozen.";
    case "WING_ONLY"
        definition.bodyCase = "wing_only";
        definition.subcase = "wing_gust_load_alleviation";
        definition.referenceType = "wing_gla";
        definition.gustEnabled = true;
        definition.executionImplemented = false;
        definition.qualificationStatus = "QUALIFIED_PREREQUISITE_ADAPTER_PENDING";
        definition.limitation = ( ...
            "The wing-only prerequisite is closed, but its public execution " + ...
            "adapter is not yet promoted through this facade.");
end
end
