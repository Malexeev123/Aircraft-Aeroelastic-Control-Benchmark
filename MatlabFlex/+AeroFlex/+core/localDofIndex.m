function idxLocal = localDofIndex(keepDofs, dofSet)
% localDofIndex
%   given e.g. dofSet = [6*(k-1)+(1:6)], find them in keepDofs 
%   returns array of local indices 
idxLocal = [];
for dd = dofSet(:).'
   iPos = find( keepDofs==dd, 1, 'first' );
   if ~isempty(iPos)
       idxLocal(end+1) = iPos; %#ok<AGROW>
   end
end