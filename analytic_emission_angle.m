function [B] = analytic_emission_angle(emission_angle, unidirectional_rms_slope, B0)
%
% This function calculates the directional emissivity of a rough Gaussain
% surface illuminated from zenith by analytically solving the intergral
% appearing in equation 41 of Smith, B., Lunar Surface Roughness: Shadowing 
% and Thermal Emission.
% For more information, see Rubanenko et al. 2020, JGR.
%
% Inputs:
% emission_angle:           the emission (observation) angle (scalar or vector)
% unidirectional_rms_slope: the unidirectional rms slope of the surface,
%                           given by the bidirectional_rms_slope / sqrt(2)
% B0:                       The infrared brightness at zero emission angle, defaults to 1 if left empty
%
% Outputs:
% B:                        The infrared brightness at some observation angle.
%
% Written by Lior Rubanenko, 2019, UCLA
%
%%%
% Example for a surface with bidirectional rms slope = 75 degrees: 
% psi = linspace(0,pi/2); urms = tand(75)/sqrt(2); 
% B_psi = analytic_emission_angle(linspace(0,pi/2), tand(75)/sqrt(2))
% plot(rad2deg(psi), B_psi); 
% xlabel('Emission angle'); ylabel(Normalized infrared brightness')
%%%

if nargin == 2
    B0 = 1
end

I_f = sqrt(pi/2./unidirectional_rms_slope.^2) .* exp(1./2./unidirectional_rms_slope.^2) ...
    .* erfc(1./sqrt(2*unidirectional_rms_slope.^2));

lambda_smith_f = @(emission_angle) 1/2 .* ((2./pi)^(1./2) * unidirectional_rms_slope./cot(emission_angle)...
    .* exp(-(cot(emission_angle)).^2./2./unidirectional_rms_slope.^2) - erfc(cot(emission_angle)./sqrt(2)./unidirectional_rms_slope));

Bpi2_f = sqrt(2*pi/unidirectional_rms_slope.^2) .* B0 ./ 4./ pi./ I_f./ unidirectional_rms_slope.^2 .* exp(1./4./unidirectional_rms_slope.^2)...
    .* (besselk(1,1./4./unidirectional_rms_slope.^2) - besselk(0,1./4./unidirectional_rms_slope.^2));

analytic_2_f = @(emission_angle) -(1+cot(emission_angle).^2) .* tan(emission_angle)...
    ./ 2 .* exp((1-cot(emission_angle).^2)./4./unidirectional_rms_slope.^2) ...
    .* (besselk(0, (1+cot(emission_angle).^2)./4./unidirectional_rms_slope.^2) - besselk(1, (1+cot(emission_angle).^2)./4./unidirectional_rms_slope.^2));

Btilde_f = @(emission_angle) B0 ./ 2./pi./I_f./unidirectional_rms_slope.^2 .* analytic_2_f(emission_angle)./(1+lambda_smith_f(emission_angle));

B = B0 - (B0 - Bpi2_f)./Bpi2_f .* Btilde_f(emission_angle);
