function [wrench, info] = recoverRootWrench(phiRoot, Pr, reactionChannel)
%RECOVERROOTWRENCH Recover the symmetric full-wing root wrench by force duality.
%   The root motion basis satisfies nuRoot = phiRoot*qdot. Virtual power
%   therefore gives reactionChannel = Pr*phiRoot.'*wrench. The full-span
%   Pazy model has range(Pr) = range(phiRoot.'). The complete dual is solved
%   before applying the physical symmetry projector. Production retains
%   [Fx,Fz,My]; lateral/roll/yaw entries are zero because asymmetric recovery
%   is not validated for the symmetry-reduced coupled plant.

arguments
    phiRoot double
    Pr double
    reactionChannel double
end

reactionChannel = reactionChannel(:);
nModal = size(phiRoot, 2);
if size(phiRoot, 1) ~= 6 || ~isequal(size(Pr), [nModal, nModal]) || ...
        numel(reactionChannel) ~= nModal
    error('AeroFlex:beam:recoverRootWrench:Dimensions', ...
        'Expected phiRoot 6-by-N, Pr N-by-N, and reactionChannel N-by-1.');
end
if any(~isfinite(phiRoot), 'all') || any(~isfinite(Pr), 'all') || ...
        any(~isfinite(reactionChannel))
    error('AeroFlex:beam:recoverRootWrench:Nonfinite', ...
        'Reaction-recovery inputs must be finite.');
end

calibration = Pr*phiRoot.';
longitudinalAxes = [1, 3, 5];
singularValues = svd(calibration);
rankTolerance = max(size(calibration))*eps(singularValues(1));
numericalRank = sum(singularValues > rankTolerance);
if numericalRank ~= 6
    error('AeroFlex:beam:recoverRootWrench:RankDeficient', ...
        'Root-wrench calibration rank is %d; six axes are required.', ...
        numericalRank);
end

wrenchFull = calibration\reactionChannel;
residual = calibration*wrenchFull - reactionChannel;
wrench = zeros(6, 1);
wrench(longitudinalAxes) = wrenchFull(longitudinalAxes);

if nargout > 1
    info = struct( ...
        'calibration', calibration, ...
        'representedAxes', longitudinalAxes, ...
        'rank', numericalRank, ...
        'singularValues', singularValues, ...
        'conditionNumber', singularValues(1)/singularValues(end), ...
        'residual', residual, ...
        'residualNorm', norm(residual));
end
end
