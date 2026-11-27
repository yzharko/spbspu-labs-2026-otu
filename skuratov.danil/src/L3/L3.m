clc;
clear;
close all;
s = tf('s');
K0 = 1.9;
T = 0.7;
T1 = 1.1;
T2 = 1.3;
T1_5 = 3;
omega = logspace(-2, 2, 1600);
K_values = [0.95 1.9 3.8];
T_values = [0.35 0.7 1.4];
output_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

names = {'Интегрирующее', 'Идеальное дифференцирующее', 'Усилительное', ...
    'Апериодическое 1-го порядка', 'Апериодическое 2-го порядка', ...
    'Колебательное', 'Консервативное', 'Инерционное интегрирующее', ...
    'Реальное дифференцирующее', 'Форсирующее', 'Изодромное'};
W = {K0/s, K0*s, tf(K0), K0/(T*s + 1), ...
    K0/(T2^2*s^2 + T1_5*s + 1), K0/(T2^2*s^2 + T1*s + 1), ...
    K0/(T*s^2 + 1), K0/(s*(T*s + 1)), K0*s/(T*s + 1), ...
    K0*(T*s + 1), K0*(T*s + 1)/s};

for number = 1:numel(W)
    models = cell(1, 3);
    labels = cell(1, 3);
    pole_values = zeros(1, 3);
    for j = 1:3
        models{j} = W{number} * K_values(j) / K0;
        labels{j} = sprintf('K = %.3g', K_values(j));
    end
    if number == 7
        pole_values(:) = 1/sqrt(T);
    end
    filename = fullfile(output_dir, sprintf('link%02d.png', number));
    draw_characteristics(models, omega, labels, names{number}, filename, number, pole_values);
end

for number = [4 7 8 9 10 11]
    pole_values = zeros(1, 3);
    for j = 1:3
        current_T = T_values(j);
        switch number
            case 4
                models{j} = K0/(current_T*s + 1);
            case 7
                models{j} = K0/(current_T*s^2 + 1);
                pole_values(j) = 1/sqrt(current_T);
            case 8
                models{j} = K0/(s*(current_T*s + 1));
            case 9
                models{j} = K0*s/(current_T*s + 1);
            case 10
                models{j} = K0*(current_T*s + 1);
            case 11
                models{j} = K0*(current_T*s + 1)/s;
        end
        labels{j} = sprintf('T = %.3g', current_T);
    end
    filename = fullfile(output_dir, sprintf('exp%02d_T.png', number));
    draw_characteristics(models, omega, labels, [names{number} ': влияние T'], filename, number, pole_values);
end

numbers = [5 5 6 6];
parameters = {'T1', 'T2', 'T1', 'T2'};
time_values = [2.6 3 4.5; 0.65 1.3 1.45; 0.55 1.1 2.2; 0.65 1.3 2.6];
for experiment = 1:4
    number = numbers(experiment);
    for j = 1:3
        current_T1 = T1;
        current_T2 = T2;
        if number == 5
            current_T1 = T1_5;
        end
        if strcmp(parameters{experiment}, 'T1')
            current_T1 = time_values(experiment, j);
        else
            current_T2 = time_values(experiment, j);
        end
        models{j} = K0/(current_T2^2*s^2 + current_T1*s + 1);
        labels{j} = sprintf('%s = %.3g', parameters{experiment}, time_values(experiment, j));
    end
    filename = fullfile(output_dir, sprintf('exp%02d_%s.png', number, parameters{experiment}));
    caption = [names{number} ': влияние ' parameters{experiment}];
    draw_characteristics(models, omega, labels, caption, filename, number, zeros(1, 3));
end

function draw_characteristics(models, omega, labels, caption, filename, number, pole_values)
    fig = figure('Visible', 'off', 'Color', 'white', 'Position', [100 100 1000 600]);
    layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax1 = nexttile(layout, 1);
    ax2 = nexttile(layout, 3);
    ax3 = nexttile(layout, 2, [2 1]);
    set([ax1 ax2], 'XScale', 'log');
    hold(ax1, 'on'); hold(ax2, 'on'); hold(ax3, 'on');
    colors = lines(3);
    gain_max = 0;
    for j = 1:3
        frequency = omega;
        valid = true(size(frequency));
        if number == 7
            frequency = sort([omega, pole_values(j)*[0.999 1 1.001]]);
            valid = frequency ~= pole_values(j);
            gain_max = max(gain_max, abs(dcgain(models{j})));
        end
        [mag, phase] = bode(models{j}, frequency(valid));
        [re, im] = nyquist(models{j}, frequency(valid));
        level = nan(size(frequency));
        phase_line = level; real_line = level; imag_line = level;
        level(valid) = 20*log10(mag(:));
        phase_line(valid) = phase(:);
        real_line(valid) = re(:);
        imag_line(valid) = im(:);
        if number == 7
            phase_line(phase_line > 90) = phase_line(phase_line > 90) - 360;
        end
        semilogx(ax1, frequency, level, 'Color', colors(j,:), 'LineWidth', 1.5);
        semilogx(ax2, frequency, phase_line, 'Color', colors(j,:), 'LineWidth', 1.5);
        if number == 3
            plot(ax3, real_line(1), imag_line(1), 'o', 'Color', colors(j,:), 'MarkerFaceColor', colors(j,:));
        else
            plot(ax3, real_line, imag_line, 'Color', colors(j,:), 'LineWidth', 1.5);
        end
    end
    xlabel(ax1, '\omega, рад/с'); ylabel(ax1, 'L, дБ'); title(ax1, 'ЛАЧХ');
    xlabel(ax2, '\omega, рад/с'); ylabel(ax2, '\phi, град'); title(ax2, 'ЛФЧХ');
    xlabel(ax3, 'Re W(j\omega)'); ylabel(ax3, 'Im W(j\omega)'); title(ax3, 'АФЧХ, \omega > 0');
    xlim(ax1, [omega(1) omega(end)]); xlim(ax2, [omega(1) omega(end)]);
    if ismember(number, [1 2])
        xlim(ax3, [-1 1]);
    elseif ismember(number, [3 10])
        xlim(ax3, [0 2*max(cellfun(@dcgain, models))]);
    end
    if number == 7
        xlim(ax3, [-10*gain_max 10*gain_max]); ylim(ax3, [-1 1]);
        title(ax3, {'АФЧХ, \omega > 0', 'Фрагмент годографа'});
    end
    grid(ax1, 'on'); grid(ax2, 'on'); grid(ax3, 'on');
    set([ax1 ax2 ax3], 'FontName', 'DejaVu Sans', 'FontSize', 16);
    key = legend(ax1, labels, 'Orientation', 'horizontal', 'FontSize', 14);
    key.Layout.Tile = 'north';
    title(layout, caption, 'Interpreter', 'none', 'FontName', 'DejaVu Sans', 'FontSize', 16);
    print(fig, filename, '-dpng', '-r180');
    close(fig);
end
