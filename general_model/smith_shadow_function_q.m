function [S] = smith_shadow_function_q(unidirectional_rms_slope_angle, solar_incidence, p, q)
%
% Calculates the probability a 1-d slope q is illuminated by a light source
% with incidence angle solar_incidence, assuming the surface slopes are
% distribution Gaussian
%
% Input:
% - rms_slope_angle: the unidirectional rms slope of the gaussian surface (deg)
% - solar_incidence: the light source incidence angle (deg)
% - q: the slope, tan(slope_angle)
%
% Output:
% S - the shadow function
%
%  Example:
%  smith_shadow_function_q(10, 25, 0.1)
%
% Written by Lior Rubanenko based on Smith 1967, JGR Vol. 72 16
% Technion, Israel Institute of Technology
% November 2023
%

heaviside = @(x) double(x > 0);
rms_slope = tand(unidirectional_rms_slope_angle);

mu = abs(cotd(solar_incidence));

Lambda_smith = @(mu) 1/2 .* (sqrt(2./pi) * rms_slope./mu ...
    .* exp(-(mu).^2./2./rms_slope.^2) - erfc(mu./sqrt(2)./rms_slope));
G_smith = @(mu) 1 ./ (1 + Lambda_smith((mu)));

if (solar_incidence > 0)
    S = heaviside(mu - p) .* G_smith(mu);
else
    S = heaviside(mu + p) .* G_smith(mu);
end
