function [cosIncidence] = cosIncidence(zenith_angle, azimuth_angle, p, q)
% The cosine of the angle between the observer and slope normal
cosIncidence = (-p .* sind(zenith_angle) .* cosd(azimuth_angle) ...
    -q .* sind(zenith_angle) .* sind(azimuth_angle) ...
    + cosd(zenith_angle)) ./ sqrt(1 + p.^2 + q.^2);

end
