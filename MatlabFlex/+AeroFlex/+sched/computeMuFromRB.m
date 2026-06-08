function mu = computeMuFromRB(rb, cfg, gust)
%COMPUTEMUFROMRB Compute scheduler coordinates from rigid-body state.
%
% Body axes convention: x forward, y right, z down.  Vertical gust is taken
% as positive body-z relative wind correction if supplied as a scalar.
if nargin < 3 || isempty(gust), gust = 0; end

if isempty(rb) || ~isfield(rb,'v_B') || norm(rb.v_B) < 1e-9
    U = cfg.flight.U_inf;
    alphaDeg = cfg.flight.aoa_deg;
    mu = [U, alphaDeg];
    return
end

vB = rb.v_B(:);
if isscalar(gust)
    vRel = [vB(1); vB(2); vB(3) - gust];
else
    g = gust(:);
    if numel(g) >= 3
        vRel = vB - g(1:3);
    else
        vRel = vB;
    end
end

U = norm(vRel);
if U < 1e-9
    U = cfg.flight.U_inf;
    alpha = deg2rad(cfg.flight.aoa_deg);
else
    alpha = atan2(vRel(3), vRel(1));
end
mu = [U, rad2deg(alpha)];
end
