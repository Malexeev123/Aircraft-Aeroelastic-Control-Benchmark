% function name = pointName(mu)
% %POINTNAME Deterministic folder-safe name for a library point.
% U = mu(1);
% a = mu(2);
% name = sprintf('pt_U%08.3f_alpha%+08.3f', U, a);
% name = strrep(name,'+','p');
% name = strrep(name,'-','m');
% name = strrep(name,'.','p');
% end
function run_id = pointName(U_inf, alpha_deg)
% function run_id = pointName(mu)
%POINTNAME Deterministic ROM-library folder name.
%
% Examples:
%   pointName(40, -1)    -> 'pt_U040_alpha_m01'
%   pointName(40,  2)    -> 'pt_U040_alpha_p02'
%   pointName(40.5,-1.5) -> 'pt_U040p5_alpha_m01p5'
    % U_inf = mu(1);
    % alpha_deg= mu(2);
    validateattributes(U_inf, {'numeric'}, {'scalar','finite','real','positive'});
    validateattributes(alpha_deg, {'numeric'}, {'scalar','finite','real'});

    Ustr = localUnsignedToken(U_inf, 3, 1);
    Astr = localSignedToken(alpha_deg, 2, 1);

    run_id = sprintf('pt_U%s_alpha_%s', Ustr, Astr);
end

function s = localUnsignedToken(x, intWidth, nDec)
    if abs(x - round(x)) < 1e-10
        fmt = sprintf('%%0%d.0f', intWidth);
    else
        fmt = sprintf('%%0%d.%df', intWidth + 1 + nDec, nDec);
    end

    s = sprintf(fmt, x);
    s = strrep(s, '.', 'p');
end

function s = localSignedToken(x, intWidth, nDec)
    if x < 0
        prefix = 'm';
    else
        prefix = 'p';
    end

    ax = abs(x);

    if abs(ax - round(ax)) < 1e-10
        fmt = sprintf('%%0%d.0f', intWidth);
    else
        fmt = sprintf('%%0%d.%df', intWidth + 1 + nDec, nDec);
    end

    s = sprintf(fmt, ax);
    s = strrep(s, '.', 'p');
    s = [prefix, s];
end