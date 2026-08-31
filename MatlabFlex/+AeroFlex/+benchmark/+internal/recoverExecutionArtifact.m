function artifact = recoverExecutionArtifact(runRoot)
%RECOVEREXECUTIONARTIFACT Load the newest complete runner MAT artifact.

arguments
    runRoot (1,1) string
end

files = dir(fullfile(runRoot,"*.mat"));
files = files(~contains(string({files.name}), ...
    ["CHECKPOINT","RAW_RUNTIME"],IgnoreCase=true));
[~,order] = sort([files.datenum],"descend");
artifact = struct("available",false,"path","", ...
    "summary",struct(),"history",struct());
for index = order
    path = string(fullfile(files(index).folder,files(index).name));
    variables = whos("-file",path);
    names = string({variables.name});
    if ~ismember("summary",names) || ~ismember("history",names)
        continue
    end
    loaded = load(path,"summary","history");
    artifact.available = true;
    artifact.path = path;
    artifact.summary = loaded.summary;
    artifact.history = loaded.history;
    return
end
end
