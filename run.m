clear; clf

% Parameters
solar_zenith_angles = [0 35 70];
observation_angle = linspace(-89.9, 89.9, 25);
observation_azimuths = [0 30 60];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Choose Scenario %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
scenario = 'rubanenko';

% Model parameters
if strcmp(scenario, 'smith')
    slope_distribution = 'gaussian';
    unid_rms_slope_angle = [55];
    hurst_exponent = 0.75; % doesn't matter in this scenario
    solar_constant = 1370;
    albedo = 0;
    emissivity = 1;
    plot_spectrum = 0;
    scattering_model = 'none';

elseif strcmp(scenario, 'rubanenko')
    slope_distribution = 'gaussian_mixture';
    unid_rms_slope_angle = [15 55];
    hurst_exponent = 0.75;
    solar_constant = 1370;
    albedo = 'foote_11';
    emissivity = 'keihm';
    plot_spectrum = 0;
    scattering_model = 'aha';
end

% Figure parameters
figure('Position',[100 100 1000 350]);
subplot_dims = [1, length(observation_azimuths)];
alpha_value = 0.15;
colormap_name = 'viridis';
colors = (colormap(colormap_name));
color_indices = round(linspace(20, size(colors,1)-20, length(solar_zenith_angles)));

% Create each panel
legend_entries = cell(length(solar_zenith_angles), 1);
tiledlayout(subplot_dims(1), subplot_dims(2), "TileSpacing","compact")
for panel = 1:length(observation_azimuths)
    nexttile
    hold on;
    box on;

    zenith_flag = true;
    for jj = 1:length(solar_zenith_angles)
        mean_radiance = zeros(size(observation_angle));
        for ii = 1:length(observation_angle)
            if observation_angle(ii) >= 0
                [mean_radiance_buff, rms_of_pdf] = radiance_gaussian_surface(slope_distribution, ...
                    unid_rms_slope_angle, ...
                    solar_zenith_angles(jj), ...
                    abs(observation_angle(ii)), observation_azimuths(panel), ...
                    'albedo', albedo, 'solar_constant', solar_constant, ...
                    'emissivity', emissivity, 'scattering_model', scattering_model, ...
                    'hurst_exponent', hurst_exponent);
            else
                [mean_radiance_buff, rms_of_pdf] = radiance_gaussian_surface(slope_distribution, ...
                    unid_rms_slope_angle, ...
                    solar_zenith_angles(jj), ...
                    abs(observation_angle(ii)), ...
                    wrapTo180(observation_azimuths(panel) + 180), ...
                    'albedo', albedo, 'solar_constant', solar_constant, ...
                    'emissivity', emissivity, 'scattering_model', scattering_model, ...
                    'hurst_exponent', hurst_exponent);
            end

            mean_radiance(ii) = mean_radiance_buff;
        end
        legend_entries{jj} = sprintf('z=%d°', round(solar_zenith_angles(jj)));
        
        cc = colors(color_indices(jj),:) + [0.05 -0.05 0];
        brightness_temp = (mean_radiance ./ 5.67e-8).^0.25;
        p(jj) = plot(observation_angle, brightness_temp, 'LineWidth', 4, 'color', cc);
        
        if zenith_flag == true
            brightness_temp_0 = brightness_temp;
            zenith_flag = false;
        end
    end
    
    % Panel formatting
    xlabel('Emission Angle');
    title(sprintf('a_o = %d°', observation_azimuths(panel)), 'FontSize', 14);
    grid on;
    xlim([-90 90]);
    ylim([175 400])
    set(gca,'fontsize',16)
    axis square
   
    xticks([-90 -45 0 45 90]);
    xtl = {'-90°', '-45°', '0°', '45°', '90°'};

    % Add legend to first panel only
    if panel == 1
        ylabel('Brightness Temperature (K)');
        legend(p,legend_entries, 'Location', 'best', 'FontSize', 16);
        xticklabels(xtl);
    elseif panel == 3
        set(gca,'YTickLabel',[])
        text(-35, 200, ['\omega=tan(',num2str(round(rms_of_pdf,1)),'\circ)'],'fontsize',16)
        xtl{1} = ''; % Hide the first (lowest) x-tick label
        xticklabels(xtl);
    else
        set(gca,'YTickLabel',[])
        xticklabels(xtl);
    end
end