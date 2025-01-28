## Generalized directional emissivity model
This model generalizes the [Smith 1967][https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/JZ072i016p04059] directional emissivity model to any illumination and observation angle, and adds a new idealized topographic model: a gaussian mixture surface.

Usage:
1. Copy all files.
2. Run run.m using Matlab (Tested on 2023b).

Examples:
```matlab
slope_distribution = 'gaussian';
unid_rms_slope_angle = [55];
hurst_exponent = 0.75; % doesn't matter in this scenario
solar_constant = 1370;
albedo = 0;
emissivity = 1;
plot_spectrum = 0;
scattering_model = 'none';
```


