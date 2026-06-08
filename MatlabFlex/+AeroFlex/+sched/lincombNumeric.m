function out = lincombNumeric(vals, w)
%LINCOMBNUMERIC Weighted sum for numeric arrays/cells of numeric arrays.
out = zeros(size(vals{1}), 'like', vals{1});
for k = 1:numel(w)
    out = out + w(k) * vals{k};
end
end
