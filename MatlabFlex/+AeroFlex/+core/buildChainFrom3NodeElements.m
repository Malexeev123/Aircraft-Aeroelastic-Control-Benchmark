function chain = buildChainFrom3NodeElements(fem)
% buildChainFrom3NodeElements
%   Reads an 8-element set of connectivities of the form:
%        e=1 => [nCur, nMid, nNext]
%        e=2 => [nNext, nX,  nY],  etc.
%   and extracts the chain of end-nodes in ascending order.

    ne = fem.num_elem;      % e.g. 8
    C  = fem.connectivities;    % [8 x 3], each row => [cur, mid, next]
    
    % Start from row 1 => first col is clamp
    node1 = C(1,1);  % e.g. clamp
    chain = node1;   % store in the chain

    for e = 1 : ne
        % row e => [nCur, nMid, nNext]
        rowE = C(e,:);
        nCur = rowE(1);
        nMid = rowE(2);
        nNxt = rowE(3);
        
        % check that nCur matches chain(end):
        if nCur ~= chain(end)
            if nCur == chain(end)+1
                disp('next wing');
                chain(end+1) = nCur;
            else
            error('Mismatch in connectivities row %d: expected nCur=%d but chain(end)=%d.', ...
                  e, nCur, chain(end));
            end
        end
        % The next node in the chain is the 3rd col => nNext
        chain(end+1:end+2) = [nNxt, nMid]; 

        % Then the next row e+1 should have first col = nNxt 
        % (which your data apparently always satisfies).
    end

    % Now chain has length = ne+1. Done.
end