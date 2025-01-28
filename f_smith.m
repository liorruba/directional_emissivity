function [f] = f_smith(unid_rms_slope_angle, solar_zenith_angle, emission_angle)
%
% Monostatic (same plane of illumination) shadowing according to
% Smith 1967 (fraction of the total projected area illuminated).
%
% Written by Lior Rubanenko based on Smith 1967, JGR Vol. 72 16
% Technion, Israel Institute of Technology
% November 2023
%

if (size(solar_zenith_angle) == 1)
    solar_zenith_angle = solar_zenith_angle * ones(size(emission_angle));
end

mu = abs(cotd(solar_zenith_angle));
mu_bar = abs(cotd(emission_angle));
rms_slope = tand(unid_rms_slope_angle); % unidirectional rms slope

heaviside = @(x) double(x > 0);

Lambda_smith = @(mu) 1/2 .* (sqrt(2./pi) * rms_slope./mu ...
    .* exp(-(mu).^2./2./rms_slope.^2) - erfc(mu./sqrt(2)./rms_slope));

G_smith = @(mu) 1 ./ (1 + Lambda_smith((mu)));

% Case 1
idx_1 = (emission_angle > solar_zenith_angle) & (emission_angle > 0);
% Case 2
idx_2 = (emission_angle <= solar_zenith_angle) & (emission_angle > 0);
% Case 3
idx_3 = (emission_angle < 0);

% Initialize f
f = zeros(size(emission_angle));

% Case 1
f(idx_1) = ones(size(solar_zenith_angle(idx_1)));
% Case 2
f(idx_2) = smith_shadow_function(unid_rms_slope_angle, solar_zenith_angle(idx_2)) .* ...
    (1 - cotd(solar_zenith_angle(idx_2))./cotd(emission_angle(idx_2))) + ...
    cotd(solar_zenith_angle(idx_2))./cotd(emission_angle(idx_2));
% Case 3
f(idx_3) = smith_shadow_function(unid_rms_slope_angle, solar_zenith_angle(idx_3))...
    .* G_smith(mu_bar(idx_3)) ...
    .* (1 + mu(idx_3) ./ mu_bar(idx_3)) ...
    + G_smith(mu(idx_3)) .* (1 - G_smith(mu_bar(idx_3))) ...
    - (mu(idx_3) ./ mu_bar(idx_3)) .* G_smith(mu_bar(idx_3));


end