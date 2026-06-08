function ROMlib = loadLibrary(pathOrStruct)
%LOADLIBRARY Load a ROM scheduler library from file or pass through a struct.
if isstruct(pathOrStruct)
    ROMlib = pathOrStruct;
    AeroFlex.sched.validateLibrary(ROMlib);
    return
end
pathOrStruct = char(pathOrStruct);
S = load(pathOrStruct);
if isfield(S,'ROMlib')
    ROMlib = S.ROMlib;
else
    error('Library file does not contain ROMlib: %s', pathOrStruct);
end
AeroFlex.sched.validateLibrary(ROMlib);
end
