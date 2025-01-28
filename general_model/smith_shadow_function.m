function [S] = smith_shadow_function(rms_slope_angle, incidence)
%
% Calculates the fraction of points illuminated on a rough Gaussian surface
% illuminated by a light source with zenith angle theta
%
% Input:
% - theta: the light source incidence angle (deg)
% - rms_slope_angle: the unidirectional rms slope of the gaussian surface (deg)
%
% Output:
% S - the shadow function
%
%  Example:
%  smith_shadow_function(10, 25)
%
% Written by Lior Rubanenko based on Smith 1967, JGR Vol. 72 16
% Technion, Israel Institute of Technology
% November 2023
%

rms_slope = tand(rms_slope_angle);

Lambda_smith = @(theta) 1/2 .* (sqrt(2./pi) * rms_slope./cotd(theta) ...
    .* exp(-(cotd(theta)).^2./2./rms_slope.^2) - erfc(cotd(theta)./sqrt(2)./rms_slope));

S = (1 - 1./2 .* erfc(cotd(incidence)./sqrt(2)./rms_slope)) ./ (Lambda_smith(incidence) + 1);
