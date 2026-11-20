clc;
clear;
close all;

s = tf('s');

K = 2.1;
T = 0.9;
T1 = 1.0;
T2 = 1.1;

T1_aperiodic = 2*T2;

W = {
    K/s
    K*s
    tf(K)
    K/(T*s + 1)
    K/(T2^2*s^2 + T1_aperiodic*s + 1)
    K/(T2^2*s^2 + T1*s + 1)
    K/(T*s^2 + 1)
    K/(s*(T*s + 1))
    K*s/(T*s + 1)
    K*(T*s + 1)
    K*(T*s + 1)/s
};

names = {
    '1. Интегрирующее'
    '2. Дифференцирующее'
    '3. Усилительное'
    '4. Апериодическое 1-го порядка'
    '5. Апериодическое 2-го порядка'
    '6. Колебательное'
    '7. Консервативное'
    '8. Инерционное интегрирующее'
    '9. Реальное дифференцирующее'
    '10. Форсирующее'
    '11. Изодромное'
};

w = logspace(-2, 2, 1000);

figure('Name', 'ЛАЧХ и ЛФЧХ всех звеньев', 'Position', [50 50 1400 900]);
layout = tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(layout, 'ЛАЧХ и ЛФЧХ типовых динамических звеньев');

for i = 1:length(W)
    [magnitude, phase] = bode(W{i}, w);
    magnitude = squeeze(magnitude);
    phase = squeeze(phase);

    ax = nexttile;

    yyaxis left;
    semilogx(w, 20*log10(magnitude), 'LineWidth', 1.2);
    ylabel('L, дБ');

    yyaxis right;
    semilogx(w, phase, '--', 'LineWidth', 1.2);
    ylabel('\phi, град');

    grid on;
    xlabel('\omega, рад/с');
    title(names{i});

    if i == 1
        legend(ax, {'ЛАЧХ', 'ЛФЧХ'}, 'Location', 'best');
    end
end

figure('Name', 'АФЧХ всех звеньев', 'Position', [50 50 1400 900]);
layout = tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(layout, 'АФЧХ типовых динамических звеньев');

for i = 1:length(W)
    [real_part, imaginary_part] = nyquist(W{i}, w);
    real_part = squeeze(real_part);
    imaginary_part = squeeze(imaginary_part);

    nexttile;

    if i == 3
        plot(real_part(1), imaginary_part(1), 'o', ...
            'MarkerFaceColor', 'auto', 'LineWidth', 1.2);
    else
        if i == 7
            resonance = 1/sqrt(T);
            resonance_index = find(w > resonance, 1);
            real_part(resonance_index) = NaN;
            imaginary_part(resonance_index) = NaN;
        end

        plot(real_part, imaginary_part, 'LineWidth', 1.2);
        hold on;
        plot(real_part, -imaginary_part, '--', 'LineWidth', 1.0);
    end

    if i == 1 || i == 2
        xlim([-1 1]);
    elseif i == 7
        xlim([-5*K 5*K]);
        ylim([-K K]);
    elseif i == 10 || i == 11
        xlim([0 2*max(abs(real_part))]);
    elseif i ~= 8
        axis equal;
    end

    grid on;
    xlabel('Re W(j\omega)');
    ylabel('Im W(j\omega)');
    title(names{i});
end

% влияние K
K_values = [0.5 2.1 5.0];
labels_K = {'K = 0.5', 'K = 2.1', 'K = 5.0'};
colors_K = lines(length(K_values));

figure('Name', 'Влияние K', 'Position', [50 50 1400 900]);
layout = tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(layout, 'Влияние коэффициента K на ЛАЧХ и ЛФЧХ всех звеньев');

for i = 1:length(W)
    ax = nexttile;
    hold on;
    set(ax, 'XScale', 'log');
    magnitude_lines = gobjects(1, length(K_values));

    yyaxis left;
    for j = 1:length(K_values)
        WK = (K_values(j)/K)*W{i};
        magnitude = bode(WK, w);
        magnitude = squeeze(magnitude);
        magnitude_lines(j) = semilogx(w, 20*log10(magnitude), ...
            'Color', colors_K(j, :), 'LineWidth', 1.2);
    end
    ylabel('L, дБ');

    % при положительном K фазовая характеристика не меняется
    [~, phase] = bode(W{i}, w);
    phase = squeeze(phase);
    yyaxis right;
    semilogx(w, phase, 'k--', 'LineWidth', 1.0);
    ylabel('\phi, град');

    grid on;
    xlabel('\omega, рад/с');
    title(names{i});

    if i == 1
        legend(ax, magnitude_lines, labels_K, 'Location', 'best');
    end
end

% влияние T
T_scale = [0.5 1 2];
labels_T = {'0.5T', 'T', '2T'};
colors_T = lines(length(T_scale));
WT = cell(length(W), length(T_scale));

for j = 1:length(T_scale)
    scale = T_scale(j);

    % в звеньях 1-3 постоянная времени отсутствует
    WT{4, j} = K/((scale*T)*s + 1);
    WT{5, j} = K/((scale*T2)^2*s^2 + (scale*T1_aperiodic)*s + 1);
    WT{6, j} = K/((scale*T2)^2*s^2 + (scale*T1)*s + 1);
    WT{7, j} = K/((scale*T)*s^2 + 1);
    WT{8, j} = K/(s*((scale*T)*s + 1));
    WT{9, j} = K*s/((scale*T)*s + 1);
    WT{10, j} = K*((scale*T)*s + 1);
    WT{11, j} = K*((scale*T)*s + 1)/s;
end

figure('Name', 'Влияние постоянных времени', 'Position', [50 50 1400 900]);
layout = tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(layout, 'Влияние постоянных времени на ЛАЧХ и ЛФЧХ всех звеньев');

for i = 1:length(W)
    ax = nexttile;

    if i <= 3
        axis off;
        text(0.5, 0.55, names{i}, 'HorizontalAlignment', 'center', ...
            'FontWeight', 'bold');
        text(0.5, 0.40, 'Параметр T отсутствует', ...
            'HorizontalAlignment', 'center');
        continue;
    end

    hold on;
    set(ax, 'XScale', 'log');
    magnitude_lines = gobjects(1, length(T_scale));

    yyaxis left;
    for j = 1:length(T_scale)
        magnitude = bode(WT{i, j}, w);
        magnitude = squeeze(magnitude);
        magnitude_lines(j) = semilogx(w, 20*log10(magnitude), '-', ...
            'Color', colors_T(j, :), 'LineWidth', 1.2);
    end
    ylabel('L, дБ');

    yyaxis right;
    for j = 1:length(T_scale)
        [~, phase] = bode(WT{i, j}, w);
        phase = squeeze(phase);
        semilogx(w, phase, '--', 'Color', colors_T(j, :), ...
            'LineWidth', 1.0);
    end
    ylabel('\phi, град');

    grid on;
    xlabel('\omega, рад/с');
    title(names{i});

    if i == 4
        legend(ax, magnitude_lines, labels_T, 'Location', 'best');
    end
end

fprintf('\nСопрягающие частоты при исходных параметрах:\n');
fprintf('Звенья 1-3: сопрягающие частоты отсутствуют.\n');
fprintf('Звено 4: omega_c = %.3f рад/с\n', 1/T);
fprintf('Звенья 5-6: omega_c = %.3f рад/с\n', 1/T2);
fprintf('Звено 7: omega_c = %.3f рад/с\n', 1/sqrt(T));
fprintf('Звенья 8-11: omega_c = %.3f рад/с\n', 1/T);

fprintf('\nСопрягающие частоты при изменении постоянных времени:\n');
for i = 1:length(T_scale)
    scale = T_scale(i);

    fprintf('%s:\n', labels_T{i});
    fprintf('  Звено 4: %.3f рад/с\n', 1/(scale*T));
    fprintf('  Звенья 5-6: %.3f рад/с\n', 1/(scale*T2));
    fprintf('  Звено 7: %.3f рад/с\n', 1/sqrt(scale*T));
    fprintf('  Звенья 8-11: %.3f рад/с\n', 1/(scale*T));
end
