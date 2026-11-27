clear;
clc;
close all;

K = 1.8;
T = 0.6;
T1 = 0.9;
T2 = 1.2;
T1_aper = 2.6; % Для 5-го звена необходимо T1_aper >= 2*T2

% Логарифмическая сетка
omega = logspace(-3, 3, 1200);

names = { ...
    'Интегрирующее (идеальное)', ...
    'Дифференцирующее (идеальное)', ...
    'Усилительное (безынерционное)', ...
    'Апериодическое 1-го порядка', ...
    'Апериодическое 2-го порядка', ...
    'Колебательное', ...
    'Консервативное', ...
    'Инерционное интегрирующее', ...
    'Реальное дифференцирующее', ...
    'Форсирующее 1-го порядка', ...
    'Изодромное'};

%% 1. Исследование частотных характеристик всех 11 типовых звеньев
W_base = buildTransferFunctions(K, T, T1, T2, T1_aper);
resp_base = {computeFrequencyResponse(W_base, omega)};

base_labels = cell(1, numel(names));
for k = 1:numel(names)
    base_labels{k} = {sprintf('K = %.2g', K)};
end

plotAllCharacteristics(omega, resp_base, base_labels, 1:11, names, 'Исходные параметры', 1.3);

%% 2. Исследование влияния изменения K при фиксированных постоянных времени
K_values = [0.5, 1.0, 1.8, 3.0, 5.0];
resp_K = cell(1, numel(K_values));
K_tile_labels = cell(1, numel(names));

for j = 1:numel(K_values)
    W_now = buildTransferFunctions(K_values(j), T, T1, T2, T1_aper);
    resp_K{j} = computeFrequencyResponse(W_now, omega);
end

K_curve_names = arrayfun(@(val) sprintf('K = %.2g', val), ...
    K_values, 'UniformOutput', false);
for k = 1:numel(names)
    K_tile_labels{k} = K_curve_names;
end

plotAllCharacteristics(omega, resp_K, K_tile_labels, 1:11, names, ...
    'Влияние коэффициента усиления K', 1.15);

%% 3. Исследование влияния изменения T при фиксированном K
% У звеньев 1-3 T отсутствует, для звеньев 4, 7-11 варьируется T,
% для звеньев 5 и 6 варьируется T1 при фиксированном T2.
T_values = T * [1/5, 1/2, 1, 2, 5];
T1_aper_values = [2.4, 2.6, 3.5, 5.0, 8.0];
T1_osc_values = [0.2, 0.5, 0.9, 1.5, 2.2];

n_T = numel(T_values);
resp_T = cell(1, n_T);
T_tile_labels = cell(1, numel(names));

for j = 1:n_T
    W_now = buildTransferFunctions(K, T_values(j), ...
        T1_osc_values(j), T2, T1_aper_values(j));
    resp_T{j} = computeFrequencyResponse(W_now, omega);
end

T_curve_names = arrayfun(@(val) sprintf('T = %.3g с', val), ...
    T_values, 'UniformOutput', false);
T1_aper_names = arrayfun(@(val) sprintf('T1 = %.2g с (\\zeta = %.2f)', ...
    val, val/(2*T2)), T1_aper_values, 'UniformOutput', false);
T1_osc_names = arrayfun(@(val) sprintf('T1 = %.2g с (\\zeta = %.2f)', ...
    val, val/(2*T2)), T1_osc_values, 'UniformOutput', false);

for k = 1:numel(names)
    if ismember(k, [1, 2, 3])
        T_tile_labels{k} = {'T отсутствует'};
    elseif k == 5
        T_tile_labels{k} = T1_aper_names;
    elseif k == 6
        T_tile_labels{k} = T1_osc_names;
    else
        T_tile_labels{k} = T_curve_names;
    end
end

resp_T_for_plot = resp_T;
for j = 2:n_T
    for k = [1, 2, 3]
        resp_T_for_plot{j}(k).L_db = nan(size(omega));
        resp_T_for_plot{j}(k).phase = nan(size(omega));
        resp_T_for_plot{j}(k).re = nan(size(omega));
        resp_T_for_plot{j}(k).im = nan(size(omega));
    end
end

plotAllCharacteristics(omega, resp_T_for_plot, T_tile_labels, 1:11, names, ...
    'Влияние постоянных времени T и T1', 1.15);

%% 4. Исследование влияния постоянной времени T2 для звеньев 5 и 6
T2_aper_values = [0.2, 0.5, 0.9, 1.2, 1.3]; % Условие T2 <= T1_aper / 2 = 1.3
T2_osc_values = [0.5, 0.8, 1.2, 2.0, 3.5];  % Условие T2 > T1 / 2 = 0.45

n_T2 = numel(T2_aper_values);
resp_T2 = cell(1, n_T2);

for j = 1:n_T2
    W_aper = buildTransferFunctions(K, T, T1, T2_aper_values(j), T1_aper);
    W_osc = buildTransferFunctions(K, T, T1, T2_osc_values(j), T1_aper);

    W_pair = {W_aper{5}, W_osc{6}};
    resp_T2{j} = computeFrequencyResponse(W_pair, omega);
end

T2_labels = { ...
    arrayfun(@(val) sprintf('T2 = %.2g с (\\zeta = %.2f)', ...
        val, T1_aper/(2*val)), T2_aper_values, 'UniformOutput', false), ...
    arrayfun(@(val) sprintf('T2 = %.2g с (\\zeta = %.2f)', ...
        val, T1/(2*val)), T2_osc_values, 'UniformOutput', false)};

plotAllCharacteristics(omega, resp_T2, T2_labels, [5, 6], names([5, 6]), ...
    'Влияние постоянной времени T2 (звенья 2-го порядка)', 1.2);

printAnalyticalSummary(K, T, T1, T2, T1_aper, K_values, T_values);

%% Вспомогателььные функции
function W = buildTransferFunctions(K, T, T1, T2, T1_aper)
    s = tf('s');

    W = cell(1, 11);
    W{1}  = K/s;
    W{2}  = K*s;
    W{3}  = tf(K);
    W{4}  = K/(T*s + 1);
    W{5}  = K/(T2^2*s^2 + T1_aper*s + 1);
    W{6}  = K/(T2^2*s^2 + T1*s + 1);
    W{7}  = K/(T*s^2 + 1);
    W{8}  = K/(s*(T*s + 1));
    W{9}  = K*s/(T*s + 1);
    W{10} = K*(T*s + 1);
    W{11} = K*(T*s + 1)/s;
end

function resp = computeFrequencyResponse(W_cell, omega)
    % bode и nyquist возвращают трёхмерные массивы, squeeze преобразует их в 1D векторы.
    n = numel(W_cell);
    resp = struct('L_db', cell(1, n), 'phase', cell(1, n), ...
                  're', cell(1, n), 'im', cell(1, n), ...
                  'isStatic', cell(1, n));

    for k = 1:n
        [mag, phase] = bode(W_cell{k}, omega);
        [re, im] = nyquist(W_cell{k}, omega);

        mag_vec = reshape(squeeze(mag), size(omega));
        phase_vec = reshape(squeeze(phase), size(omega));
        re_vec = reshape(squeeze(re), size(omega));
        im_vec = reshape(squeeze(im), size(omega));

        resp(k).L_db = 20 * log10(mag_vec);
        resp(k).phase = phase_vec;
        resp(k).re = re_vec;
        resp(k).im = im_vec;
        resp(k).isStatic = isstatic(W_cell{k});
    end
end

function plotAllCharacteristics(omega, respSets, tileLabels, ...
    blockNumbers, blockNames, groupTitle, lineWidth)

    plotFrequencyGrid(omega, respSets, tileLabels, blockNumbers, blockNames, ...
        ['ЛАЧХ: ', groupTitle], ['ЛАЧХ — ', groupTitle], 'lach', lineWidth);

    plotFrequencyGrid(omega, respSets, tileLabels, blockNumbers, blockNames, ...
        ['ЛФЧХ: ', groupTitle], ['ЛФЧХ — ', groupTitle], 'lfch', lineWidth);

    plotFrequencyGrid(omega, respSets, tileLabels, blockNumbers, blockNames, ...
        ['АФЧХ: ', groupTitle], ['АФЧХ — ', groupTitle], 'afch', lineWidth);
end

function plotFrequencyGrid(omega, respSets, tileLabels, ...
    blockNumbers, blockNames, figureName, plotTitle, charType, lineWidth)

    figure('Name', figureName, 'Color', 'w');

    n_blocks = numel(blockNumbers);
    if n_blocks <= 2
        tiledlayout(1, n_blocks, 'TileSpacing', 'compact', 'Padding', 'compact');
    else
        tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    end

    for k = 1:n_blocks
        nexttile;
        hold on;
        grid on;

        labels_k = tileLabels{k};
        n_active = numel(labels_k);

        for j = 1:n_active
            data = respSets{j}(k);
            label_now = labels_k{j};

            switch charType
                case 'lach'
                    semilogx(omega, data.L_db, 'LineWidth', lineWidth, ...
                        'DisplayName', label_now);
                    set(gca, 'XScale', 'log');
                    xlabel('Частота \omega, рад/с');
                    ylabel('L(\omega), дБ');

                case 'lfch'
                    semilogx(omega, data.phase, 'LineWidth', lineWidth, ...
                        'DisplayName', label_now);
                    set(gca, 'XScale', 'log');
                    xlabel('Частота \omega, рад/с');
                    ylabel('\phi(\omega), град');

                case 'afch'
                    if data.isStatic
                        plot(data.re(1), data.im(1), 'o', ...
                            'LineWidth', lineWidth, 'DisplayName', label_now);
                    else
                        % Сплошная ветвь соответствует w > 0, пунктирная w < 0
                        h = plot(data.re, data.im, 'LineWidth', lineWidth, ...
                            'DisplayName', label_now);
                        plot(flip(data.re), -flip(data.im), '--', ...
                            'LineWidth', lineWidth * 0.85, 'Color', h.Color, ...
                            'HandleVisibility', 'off');
                    end
                    xlabel('Re W(j\omega)');
                    ylabel('Im W(j\omega)');

                otherwise
                    error('Неизвестный тип характеристики: %s', charType);
            end
        end

        title(sprintf('%d. %s', blockNumbers(k), blockNames{k}), ...
            'Interpreter', 'none');
        legend('show', 'Location', 'best');
        hold off;
    end

    sgtitle(plotTitle);
end

% Расчёт чисел для вывода
function printAnalyticalSummary(K0, T0, T1, T2, T1_aper, K_values, T_values)
    omega_c0 = 1 / T0;
    disc_aper = sqrt(T1_aper^2 - 4*T2^2);
    Ta = (T1_aper + disc_aper) / 2;
    Tb = (T1_aper - disc_aper) / 2;
    omega_c1 = 1 / Ta;
    omega_c2 = 1 / Tb;

    omega_n6 = 1 / T2;
    zeta6 = T1 / (2 * T2);
    omega_r6 = omega_n6 * sqrt(1 - 2*zeta6^2);
    L_max6 = 20 * log10(K0 / (2 * zeta6 * sqrt(1 - zeta6^2)));
    omega_n7 = 1 / sqrt(T0);

    fprintf('Звенья 4 и 8-11: omega_c = 1/T0 = %.4f рад/с\n', omega_c0);
    fprintf('Звено 5: omega_c1 = %.4f рад/с, omega_c2 = %.4f рад/с\n', ...
        omega_c1, omega_c2);
    fprintf('Звено 6: omega_n = %.4f рад/с, zeta = %.4f, omega_r = %.4f рад/с, L_max = %.2f дБ\n', ...
        omega_n6, zeta6, omega_r6, L_max6);
    fprintf('Звено 7: omega_n = 1/sqrt(T0) = %.4f рад/с\n\n', ...
        omega_n7);

    % Влияние изменения T на сопрягающие частоты
    fprintf('Влияние постоянной времени T на сопрягающую частоту 4, 8-11\n');
    fprintf('   T (с)     omega_c = 1/T рад/с\n');
    for Tv = T_values
        fprintf('%8.3f          %10.4f\n', Tv, 1/Tv);
    end

    % Влияние изменения K на сдвиг ЛАЧХ и точки пересечения осей
    fprintf('\nВлияние коэффициента K на сдвиг ЛАЧХ и частоты среза\n');
    fprintf('   K       Delta L   L(0)  w_ср зв.1 рад/с   w_ср зв.4 рад/с\n');
    for Kv = K_values
        delta_L = 20 * log10(Kv / K0);
        L0_db = 20 * log10(Kv);
        w_cross_1 = Kv;
        if Kv > 1
            w_cross_4 = sqrt(Kv^2 - 1) / T0;
            fprintf('%7.2f %9.2f %8.2f %10.4f    %10.4f\n', ...
                Kv, delta_L, L0_db, w_cross_1, w_cross_4);
        else
            fprintf('%7.2f %9.2f %8.2f %10.4f           нет\n', ...
                Kv, delta_L, L0_db, w_cross_1);
        end
    end
end
