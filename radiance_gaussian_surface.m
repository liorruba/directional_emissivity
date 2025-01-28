function [mean_radiance, rms_of_pdf] = ...
    radiance_gaussian_surface (...
    slope_distribution, rms_slope_angle, ...
    solar_zenith_angle, observation_angle, observation_azimuth, ...
    varargin ...
    )
% ----
%
% The thermal radiance ("brightness") of a rough Gaussian surface. This
% model extends the model of Smith 1967 for any illumination and
% observation angle.
%
% References:
% (*) Smith, 1967. JGR, 72(16), pp.4059-4067.
% (*) Rubanenko et al., 2020. JGR, 125(6), p.e2020JE006377.
% (*) See Rubanenko & Powell, 2025 for complete model derivation.
%
% ----
%
%%%%%%%%%%
% Input: %
%%%%%%%%%%
% Required:
% - slope_distribution: the slope distribution of the surface. Can be
%   either 'gaussian' or 'gaussian_mixture'. Gaussian is a single-lobe
%   gaussian surface, as in Smith 1967. Gaussian mixture is a mixture of
%   Gaussian PDFs with difference rms slopes, as an approximation
%   for a fractal surface.
% - rms_slope_angle: the root mean square slope angle of the surface, in
%   degrees.
%   - If slope_distribution is 'gaussian': rms_slope_angle is a scalar, as
%     in Smith 1967.
%   - If slope_distribution is 'gaussian_mixture': pass the maximum and 
%     minimum rms slopes of the surface as a two-element vector [min,max]. 
%     To control roughness scaling, also pass hurst_exponent (see below).
% - solar_zenith_angle: the angle between the solar vector and the mean
%   plane, in degrees. Range: [0, 90). Default: 0 deg.
% - observation_angle: the angle between the observer and the mean plane,
%   in degrees. Range: [0, 90). Note: sometimes termed "emission angle".
%   Default: 60 deg.
% - observation_azimuth: the azimuth difference between the Sun and the
%   observer, in degrees. Range: [0, 180]. Default: 0 deg.
%
% Optional:
% - solar_constant: the mean solar flux at the planetary surface. Default:
%   1370 W m^-2.
% - hurst_exponent: only works when slope_distribution is 'dualscale_gaussian'.
%   Sets scaling between the roughness magnitude at different lateral
%   scales.
% - scattering_model: emitted radiation from topography.
%   - If 'none', doesn't include scattering.
%   - If 'aha', uses the model of Aharonson & Schorghofer 2006.
%   - If 'buhl', uses the model of Buhl 1968.
% - albedo: the surface albedo (reflectance for solar radiation).
%   - If scalar: constant scalar value.
%   - If 'keihm':
%     uses the model of Keihm (1984), with parameters by Keihm (1984):
%     A = A0 + a .* (theta./45).^3 + b .* (theta./90).^8, where theta is
%     the solar incidence angle and where A0 = 0.12; a = 0.06; b = 0.25.
%   - If 'hayne': uses the Keihm (1984) model with parameters by Hayne et al.,
%     2017, with A0 = 0.12; a = 0.03; b = 0.14.
%   - If 'foote_11', uses the Foote et al. (2020) albedo model:
%     A = A0 + a .* theta.^2 + b .* theta, with parameters fitted using
%     data from Apollo 11.
%     A0 = 0.068061; a = 9.4032e-6; b = -1.8345e-4;
%   - If 'foote_16', uses the Foote et al. (2020) albedo model, with
%     parameters fitted using Apollo 16 data:
%     A0 = 1.44e-1; a = 1.33e-5; b = -1.025e-4;
%     Default: 0.
% - emissivity: the surface thermal emissivity.
%   - If scalar, uses a constant emissivity.
%   - If 'keihm', uses the Keihm (1984) temperature dependent emissivity
%     model
%     Default: 1.
%
%%%%%%%%%%%
% Output: %
%%%%%%%%%%%
% - mean_radiance: the mean radiance of the rough surface (W / m^2 / sr).
%
%%%%%%%%%%%
% Output: %
%%%%%%%%%%%
%
% Written by Lior Rubanenko
% Technion, Israel Institute of Technology
% Planetary Science Intitute
% November 2023
% Updated Decemeber 2024
%
%%%%
% Cite: Rubanenko & Powell, 2025.
%%%%

%% Check input using the inputparser object:
% Initialize the inputparser object:
p = inputParser;

% Set default values:
def.hurst_exponent = 0.5;
def.scattering_model = 'none';
def.albedo = 0;
def.emissivity = 1;
def.solar_constant = 1370;
def.plot_spectrum = 'false';

% Add input verification:
addRequired(p, 'slope_distribution', @(x) ischar(x) || isstring(x));  % Allow both char and string
addRequired(p, 'rms_slope_angle', @isnumeric);
addRequired(p, 'solar_zenith_angle', @isnumeric);
addRequired(p, 'observation_angle', @isnumeric);
addRequired(p, 'observation_azimuth', @isnumeric);

% Add optional parameters with their default values
addParameter(p, 'hurst_exponent', def.hurst_exponent, @isnumeric);
addParameter(p, 'scattering_model', def.scattering_model, @(x) ischar(x) || isstring(x));
addParameter(p, 'albedo', def.albedo, @(x) isnumeric(x) || ischar(x) || isstring(x));
addParameter(p, 'emissivity', def.emissivity, @(x) isnumeric(x) || ischar(x) || isstring(x));
addParameter(p, 'plot_spectrum', def.plot_spectrum, @islogical);
addParameter(p, 'solar_constant', def.solar_constant, @isnumeric);

% Parse input arguments:
parse(p, slope_distribution, rms_slope_angle, ...
    solar_zenith_angle, observation_angle, observation_azimuth, ...
    varargin{:});  % Added {:} to properly expand varargin

% Extract results:
hurst_exponent = p.Results.hurst_exponent;
scattering_model = p.Results.scattering_model;
albedo = p.Results.albedo;
emissivity = p.Results.emissivity;
solar_constant = p.Results.solar_constant;
plot_spectrum = p.Results.plot_spectrum;

%%%%%%%%%%%%%
% Constants %
%%%%%%%%%%%%%
sb_const = 5.67e-8; % Stefan Boltzmann constant in SI

%%%%%%%%%%%%%%%%%
% Albedo models %
%%%%%%%%%%%%%%%%%
if strcmp(albedo,  "powell")
    A0 = 0.14; a = 0.007; b = 0.0152;
    albedo_powell = @(theta) A0 + a .* (theta./45).^3 + b .* (theta./90).^8;
    albedo = @(p,q) albedo_powell(acosd(cosIncidence(observation_angle, observation_azimuth, p, q)));

elseif strcmp(albedo,  "keihm")
    A0 = 0.12; a = 0.03; b = 0.14;
    albedo_keihm = @(theta) A0 + a .* (theta./45).^3 + b .* (theta./90).^8;
    albedo = @(p,q) albedo_keihm(acosd(cosIncidence(observation_angle, observation_azimuth, p, q)));

elseif strcmp(albedo,  "hayne")
    A0 = 0.12; a = 0.06; b = 0.25;
    albedo_keihm = @(theta) A0 + a .* (theta./45).^3 + b .* (theta./90).^8;
    albedo = @(p,q) albedo_keihm(acosd(cosIncidence(observation_angle, observation_azimuth, p, q)));

elseif strcmp(albedo,  "foote_16") % Define inc. angle dependent albedo (Foote & Paige 2020):
    A0 = 1.44e-1; a = 1.33e-5; b = -1.025e-4;
    albedo_foote = @(theta) A0 + a .* theta.^2 + b .* theta;
    albedo = @(p,q) albedo_foote(acosd(cosIncidence(observation_angle, observation_azimuth, p, q)));

elseif strcmp(albedo,  "foote_11") % Define inc. angle dependent albedo (Foote & Paige 2020):
    A0 = 0.068061; a = 9.4032e-6; b = -1.8345e-4;
    albedo_foote = @(theta) A0 + a .* theta.^2 + b .* theta;
    albedo = @(p,q) albedo_foote(acosd(cosIncidence(observation_angle, observation_azimuth, p, q)));

elseif isscalar(albedo)
    albedo = @(p,q) albedo;

else
    error('ERROR: albedo can be one of ''powell'', ''keihm'', ''hayne'', ''foote_11/16'' or a scalar.')
end

%%%%%%%%%%%%%%%%%%%%%
% Emissivity models %
%%%%%%%%%%%%%%%%%%%%%
if strcmp(emissivity,  "keihm") % Uses the emissivity measured by Birkebak 1970 (Keihm 1984)
    eps0 = 0.9696;
    eps1 = 0.9664e-4;
    eps2 = -0.3167e-6;
    eps3 = -0.50691e-9;

    emissivity_birk = @(T) eps0 + eps1 .* T + eps2 .* T.^2 + eps3 .* T.^3;
    emissivity = @(T) emissivity_birk(T);

elseif isscalar(emissivity)
    emissivity = @(T) emissivity;

else
    error('ERROR: emissivity can be ''keihm'' or a scalar.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slope distribution models %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The rms slope of the surface is the tangent of the rms slope angle:
rms_slope = tand(rms_slope_angle);

if strcmp(slope_distribution, 'gaussian')
    % The Gaussian probability density function for the slope vector components
    % p and q:
    if length(rms_slope) ~= 1
        error("ERROR: when slope_distribution is 'gaussian', rms_slope must be a scalar.")
    end

    PDF = @(p,q) exp(-(p.^2 + q.^2) / 2 / rms_slope.^2) / 2 / pi / rms_slope.^2;
    rms_of_pdf = atand(rms_slope);
    
    % t_vis_illum is defined in Smith 1967 as T: the probability of a facet to
    % be seen by both the observer and the Sun.
    t_vis_illum = bistatic_visibility_function(rms_slope_angle, solar_zenith_angle, ...
        observation_angle, observation_azimuth);

elseif strcmp(slope_distribution, 'gaussian_mixture')
    if length(rms_slope) ~= 2
        
    end
    % The Gaussian probability density function for the slope vector components
    % p and q:
    [PDF, total_rms_slope, rms_slopes, weights] = gaussian_mixture_rough_surface(rms_slope(1), rms_slope(2), hurst_exponent, 10);
    
    rms_of_pdf = atand(total_rms_slope ./ sqrt(2));
        
    % t_vis_illum is defined in Smith 1967 as T: the probability of a facet to
    % be seen by both the observer and the Sun.   
    t_vis_illum = bistatic_visibility_function(rms_slope_angle, solar_zenith_angle, ...
            observation_angle, observation_azimuth);

    t_vis_illum = @(p,q) 0;
    for rr=1:length(rms_slopes)
        % The
        f = bistatic_visibility_function(atand(rms_slopes(rr).*sqrt(2)), solar_zenith_angle, ...
            observation_angle, observation_azimuth);
        t_vis_illum = @(p,q) t_vis_illum(p,q) + weights(rr) .* f(p,q);
    end

else

    error(['ERROR: slope_distribution can be one of ''gaussian''', ...
        '''gaussian_mixture'''.'])
end

%%%%%%%%%%%%%%%%%%%%%
% Scattering models %
%%%%%%%%%%%%%%%%%%%%%
if strcmp(scattering_model, 'none')
    scattered_flux = @(p,q) 0;
elseif strcmp(scattering_model, 'aha') % Scattering per Aharonson & Schorghofer, 2006
    scattered_flux = @(p,q) scattering_aha(p, q, solar_constant, solar_zenith_angle, emissivity, albedo(0,0));
elseif strcmp(scattering_model, 'buhl')
    scattered_flux = @(p,q) scattering_buhl(p, q, rms_slope, solar_constant, solar_zenith_angle, emissivity, albedo(0,0));
else
    error('ERROR: supported scattering models are ''aha'' or ''buhl''');
end

%%%%%%%%%%%%%%%%%%
% Radiance model %
%%%%%%%%%%%%%%%%%%
% Solar flux
solar_flux = @(p,q) solar_constant ...
    .* (1 - albedo(p,q)) ...
    .* cosIncidence(solar_zenith_angle, 0, p, q);

total_flux = @(p,q) solar_flux(p,q) .* t_vis_illum(p,q) + ...
    scattered_flux(p,q);

% Emissivity correction:
% Calculate temperatures with emissivity = 1:
temp_eps1 = @(p,q) (total_flux(p,q) ./ sb_const).^0.25;

% Correct the temperatures. This computes T^4 to use later:
emissivity_corrected_temperatures_4 = @(p,q) (emissivity(temp_eps1(p, q)) .* temp_eps1(p, q).^4);

% Radiance of a facet with slope (p,q) (radiative Eq.)
L_pq = @(p, q) sb_const .* emissivity_corrected_temperatures_4(p,q) ...
    .* cosIncidence(observation_angle, observation_azimuth, p, q) ...
    ./ cosd(observation_angle) ...
    ./ cos_slope(p, q);

% The area of the surface, projected on the mean plane
I = integral2(@(p,q) PDF(p, q) .* cos_slope(p, q), -inf, inf, -inf, inf);

% Compute the mean radiance integral
mean_radiance = integral2(@(p,q) PDF(p, q) .* L_pq(p, q) , -inf, inf, -inf, inf) ./ I;



%% Helper functions
% Scattering using Aharonson & Schroghofer 2006 model
    function scattered_flux = scattering_aha(p, q, solar_constant, solar_zenith_angle, emissivity, albedo)
        sigma_T_land_4 = (solar_constant .* cosd(solar_zenith_angle));

        % Include emissivity to compute the scattered flux:
        eps_sigma_T_land_4 = emissivity(1) .* sigma_T_land_4;

        scattered_flux = eps_sigma_T_land_4 .* (1 - cos_slope(p,q)) / 2;
    end

% Scattering using Buhl 1968
    function scattered_flux = scattering_buhl(p, q, rms_slope, solar_constant, solar_zenith_angle, emissivity, albedo)
        depth_to_diameter = 2 .* (1 - cos(rms_slope));
        % depth_to_diameter = rms_slope/2;
        scattered_T = sphericalCraterTemperature(solar_constant, depth_to_diameter, albedo, emissivity(1), 90 - solar_zenith_angle);

        scattered_flux = emissivity(1) .* sb_const .* scattered_T.^4;
    end


% The surface slope angle, expressed using the slope components p and q:
    function cosslope = cos_slope(p, q)
        cosslope = 1 ./ sqrt(1 + p.^2 + q.^2);
    end

end
