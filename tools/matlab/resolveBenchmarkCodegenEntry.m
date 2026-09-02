function entryPath = resolveBenchmarkCodegenEntry(functionName,project)
%RESOLVEBENCHMARKCODEGENENTRY Resolve a packaged M-file for MATLAB Coder.
%   MATLAB R2023b requires the physical M-file path when a code-generation
%   entry point belongs to a package. Later releases also accept that form,
%   so the native builders use it consistently across supported releases.

arguments
    functionName (1,1) string
    project (1,1) struct
end

entryPath = string(which(char(functionName)));
assert(strlength(entryPath)>0 && isfile(entryPath), ...
    "AeroFlex:CodegenEntryResolution", ...
    "Cannot resolve the MATLAB Coder entry point: %s",functionName);
assert(endsWith(entryPath,".m","IgnoreCase",true), ...
    "AeroFlex:CodegenEntryType", ...
    "The MATLAB Coder entry point must resolve to an M-file: %s",entryPath);

repositoryRoot = string(project.repositoryRoot);
entryCanonical = string(java.io.File(char(entryPath)).getCanonicalPath());
rootCanonical = string(java.io.File(char(repositoryRoot)).getCanonicalPath());
assert(startsWith(entryCanonical,rootCanonical + filesep,"IgnoreCase",true), ...
    "AeroFlex:CodegenEntryScope", ...
    "The MATLAB Coder entry point is outside this repository: %s",entryPath);
end
