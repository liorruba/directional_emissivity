function [T] = bistatic_visibility_function(unid_rms_slope_angle, solar_incidence,...
    detector_incidence, detector_azimuth)
%
% The bistatic visibility (originally shadowing) function, in the case of 
% (1) same plane of illumination (from Smith 1967) or (2) a different plane of 
% illumination (Bass & Fuchs, 1975).
%
% Input:
% - rms_slope_angle: the unidirectional rms slope angle of the gaussian surface (deg)
% - solar_incidence: the solar incidence angle. Range = [0, 90]
% - detector_incidence: the incidence angle of the detector (emission
%   angle). Range = [0, 90]
% - detector_azimuth: the azimuth of the detector (emission azimuth).
%   Range = [-180, 180]
%
% Output:
% - T: the bistatic visibility function (keeping the same notation as Smith
%      1967), or the probability a slope is visibly illuminated as a function
%      of the local slopes (p,q).
%
% Written by Lior Rubanenko
% Technion, Israel Institute of Technology
% November 2023
% Updated November 2024
%

% Convert the input observer angles to conform with Smith 1967 notation
if (detector_azimuth < -180) || (detector_azimuth > 180)
    error('ERROR: detector azimuth must be in [-180, 180]')
end

if (detector_incidence < 0) || (detector_incidence >= 90)
    error('ERROR: detector incidence must be in [0, 90)')
end

if (solar_incidence < 0) || (solar_incidence >= 90)
    error('ERROR: solar incidence must be in [0, 90')
end

% Compute the rms slope from the rms slope angle
rms_slope = tand(unid_rms_slope_angle);

% Define v1, v2 as in Bourlier, 2002
v1 = abs(cotd(solar_incidence)) ./ sqrt(2) ./ rms_slope;
v2 = abs(cotd(detector_incidence)) ./ sqrt(2) ./ rms_slope;

LAMBDA = @(v) (exp(-v.^2) ./ v ./ sqrt(pi) - erfc(v)) ./ 2;

if (detector_azimuth == 0 || detector_azimuth == 180) % Single plane bistatic
    if (detector_azimuth == 180)
        detector_incidence = -detector_incidence;
    end

    if (detector_incidence < solar_incidence) && (detector_incidence > 0)
        T = @(p, q) smith_shadow_function_q(unid_rms_slope_angle, solar_incidence, p, q);

    elseif (detector_incidence >= solar_incidence)
        T = @(p, q) smith_shadow_function_q(unid_rms_slope_angle, detector_incidence, p, q);

    else
        T = @(p, q) smith_shadow_function_q(unid_rms_slope_angle, solar_incidence, p, q)...
            .* smith_shadow_function_q(unid_rms_slope_angle, detector_incidence, p, q);
    end
else % Dual plane bistatic
        T = @(p,q) P_self_shadow(solar_incidence, 0, p, q) ...
            .* P_self_shadow(detector_incidence, detector_azimuth, p, q) ...
            ./ (1 + LAMBDA(v1) + LAMBDA(v2));
end


%% Helper functions
function [P] = P_self_shadow(incidence, azimuth, p, q)
% Returns 1 if self shadowed
heaviside = @(x) double(x > 0); % Heaviside step function

P = heaviside(cosIncidence(incidence, azimuth, p, q));
end

end