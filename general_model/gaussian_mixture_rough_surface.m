function [pdfHandle, bidirectional_rms_slope, rms_slopes, weights] = ...
    gaussian_mixture_rough_surface (...
    rms_slope_min, rms_slope_max, ...
    scale_exponent, num_scales...
    )
% This function returns the probability density function (PDF) of a
% mixture of Gaussian surfaces with directional slopes p, q.
%%%%%%%%%%
% Input: %
%%%%%%%%%%
% - rms_slope_min, rms_slope_max: the minimum and maximum rms slopes of the
%   Gaussian distributions to be mixed. 
% - scale_exponent: controls the scaling of roughness. This is similar (but
%   not the same) as the hurst exponent for a true fractal surface, which
%   controls the roughness magnitude at different lateral scales.
% - num_scales: number of scales to be summed.
%
%%%%%%%%%%%
% Output: %
%%%%%%%%%%%
% - pdfHandle: a function handle f(p,q) to the PDF.
% - bidirectional_rms_slope: the bidirectional rms slope of the surface.
% - rms_slopes: individual rms slopes of each distribution.
% - weights: the weights of the mixture.
%
% Written by Lior Rubanenko
% Technion, Israel Institute of Technology
% Planetary Science Intitute
% Decemeber 2024
%
    rms_slopes = logspace(log10(rms_slope_min), log10(rms_slope_max), num_scales); % RMS slopes
    logScales = log(rms_slopes); % Log-spacing for scale dependence
    logDiffs = diff(logScales);

    intervals = zeros(1, num_scales);
    intervals(1) = logDiffs(1) / 2;
    intervals(end) = logDiffs(end) / 2;
    for i = 2:(num_scales - 1)
        intervals(i) = (logDiffs(i - 1) + logDiffs(i)) / 2;
    end

    weights = intervals .* (rms_slopes / rms_slope_min).^(-2 * scale_exponent);
    weights = weights / sum(weights); % Normalize weights

    pdfHandle = @(p, q) sumOfGaussians(p, q, rms_slopes, weights);
    bidirectional_rms_slope = sqrt(integral2(@(p,q) (p.^2+q.^2) .* pdfHandle(p,q), -inf, inf, -inf, inf));
end

%% Helper function
function val = sumOfGaussians(p, q, rmsSlopes, weights)
% This function sums the pdfs of single lobe gaussians 
    val = 0;
    for i = 1:length(rmsSlopes)
        omega = rmsSlopes(i);
        w = weights(i);
        val = val + w .* (1 / (2 * pi * omega^2)) .* exp(-(p.^2 + q.^2) / (2 * omega^2));
    end
end
