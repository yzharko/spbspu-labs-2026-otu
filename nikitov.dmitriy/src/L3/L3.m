clc;
clear;
close all;

K = 3.5;
T = 1.2;
T1 = 1.4;
T2 = 1.8;
T1a = 4.0;
s = tf('s');
W = {K/s, K*s, tf(K), K/(T*s + 1), ...
    K/(T2^2*s^2 + T1a*s + 1), K/(T2^2*s^2 + T1*s + 1), ...
    K/(T*s^2 + 1), K/(s*(T*s + 1)), K*s/(T*s + 1), ...
    K*(T*s + 1), K*(T*s + 1)/s};
w = logspace(-2, 2, 1600)';
gainFactors = [0.5, 1, 2];
gainShifts = 20*log10(gainFactors);

%% Частотные характеристики
amplitude = panels('ЛАЧХ типовых звеньев');
phase = panels('ЛФЧХ типовых звеньев');
nyquistAxes = panels('АФЧХ типовых звеньев');
for i = 1:numel(W)
    pole = 0;
    if i == 7
        pole = 1/sqrt(T);
    end
    [L, phi, valid] = bode_data(W{i}, w, pole);
    name = sprintf('№%d', i);
    log_plot(amplitude, w, L + gainShifts, name, 'L, дБ');
    yline(0, ':', 'HandleVisibility', 'off');
    if i == 1
        key = legend(compose('K = %g', K*gainFactors), 'Orientation', 'horizontal');
        key.Layout.Tile = 'south';
    end
    log_plot(phase, w, phi, name, '\phi, °');

    [re, im] = nyquist(W{i}, w(valid));
    U = nan(size(w));
    V = U;
    U(valid) = re(:);
    V(valid) = im(:);
    if i == 1 || i == 2
        U(:) = 0;
    elseif i == 7
        V(~isnan(V)) = 0;
    end
    nexttile(nyquistAxes);
    plot(U, V, 'b-', U, -V, 'r--');
    if i == 3
        hold on;
        plot(K, 0, 'bo', 'MarkerFaceColor', 'b');
    end
    if ismember(i, [1, 2, 8, 10, 11])
        ylim([-3*K, 3*K]);
    elseif i == 7
        xlim([-4*K, 4*K]);
    elseif ismember(i, [4, 9])
        axis equal;
    end
    if i == 10 || i == 11
        xlim([0, 2*max(U)]);
    end
    xlabel('Re W');
    ylabel('Im W');
    grid on;
    title(name);
end

%% Влияние постоянных времени
Tvalues = T * [0.5, 1, 2];
studies = {
    4, 'T',   Tvalues,        @(p) K/(p*s + 1)
    5, 'T_{1a}', [3.7, T1a, 6], @(p) K/(T2^2*s^2 + p*s + 1)
    5, 'T_2', [0.9, T2, 1.95], @(p) K/(p^2*s^2 + T1a*s + 1)
    6, 'T_1', T1*[0.5, 1, 2], @(p) K/(T2^2*s^2 + p*s + 1)
    6, 'T_2', T2*[0.5, 1, 2], @(p) K/(p^2*s^2 + T1*s + 1)
    7, 'T',   Tvalues,        @(p) K/(p*s^2 + 1)
    8, 'T',   Tvalues,        @(p) K/(s*(p*s + 1))
    9, 'T',   Tvalues,        @(p) K*s/(p*s + 1)
   10, 'T',   Tvalues,        @(p) K*(p*s + 1)
   11, 'T',   Tvalues,        @(p) K*(p*s + 1)/s
};
timeAmplitude = panels('Влияние постоянных времени на ЛАЧХ');
timePhase = panels('Влияние постоянных времени на ЛФЧХ');
for i = 1:size(studies, 1)
    [id, parameter, values, model] = studies{i,:};
    for j = 1:numel(values)
        pole = 0;
        if id == 7
            pole = 1/sqrt(values(j));
        end
        [familyL(:,j), familyPhase(:,j)] = bode_data(model(values(j)), w, pole);
    end
    name = sprintf('№%d, %s', id, parameter);
    log_plot(timeAmplitude, w, familyL, name, 'L, дБ');
    legend(compose('%g', values), 'Location', 'best');
    phaseLines = log_plot(timePhase, w, familyPhase, name, '\phi, °');
    if id == 7
        for j = 1:numel(values)
            xline(1/sqrt(values(j)), ':', 'Color', phaseLines(j).Color, ...
                'HandleVisibility', 'off');
        end
        xlim([0.35, 2.5]);
        ylim([-195, 15]);
    end
    legend(compose('%g', values), 'Location', 'best');
end

function tiles = panels(name)
    fig = figure('Name', name, 'NumberTitle', 'off', 'Theme', 'light', ...
        'Position', [60, 60, 850, 950]);
    tiles = tiledlayout(fig, 4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
end

function curves = log_plot(tiles, w, y, name, unit)
    nexttile(tiles);
    curves = semilogx(w, y, 'LineWidth', 1.2);
    xlabel('\omega, рад/с');
    ylabel(unit);
    title(name);
    grid on;
end

function [L, phi, valid] = bode_data(sys, w, pole)
    valid = abs(w - pole) > 0.01*pole;
    L = nan(size(w));
    phi = L;
    [mag, phaseDeg] = bode(sys, w(valid));
    L(valid) = 20*log10(mag(:));
    phi(valid) = phaseDeg(:);
    if pole > 0
        phi(valid) = -180*(w(valid) > pole);
    end
end
