function name = pointName(mu)
%POINTNAME Deterministic folder-safe name for a library point.
U = mu(1);
a = mu(2);
name = sprintf('pt_U%08.3f_alpha%+08.3f', U, a);
name = strrep(name,'+','p');
name = strrep(name,'-','m');
name = strrep(name,'.','p');
end
