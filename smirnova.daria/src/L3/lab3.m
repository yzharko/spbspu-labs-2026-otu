clc;
clear;
close all;
s = tf('s');
K = 2.4;
T = 0.4;
T1 = 0.8;
T2 = 0.9;
T1_5 = 2.0; % Минимально подходящая T1к пункту №5, из-за комплексности корней 
omega = logspace(-2, 3, 2000);

system_links = {
K / s,                                  '1. Интегрирующее';
    K * s,                              '2. Дифференцирующее';
    tf(K, 1),                           '3. Усилительное';
    K / (T * s + 1),                    '4. Апериодическое 1-го';
    K / (T2^2 * s^2 + T1_5 * s + 1),    '5. Апериодическое 2-го';
    K / (T2^2 * s^2 + T1 * s + 1),      '6. Колебательное';
    K / (T^2 * s^2 + 1),                '7. Консервативное';
    K / (s * (T * s + 1)),              '8. Инерционное инт.';
    (K * s) / (T * s + 1),              '9. Реальное диф.';
    K * (T * s + 1),                    '10. Форсирующее';
    (K * (T * s + 1)) / s,              '11. Изодромное'
};

num_links = size(system_links, 1);
figure('Name', 'ЛАЧХ и ЛФЧХ всех звеньев', 'NumberTitle', 'off');
for i = 1:num_links
    subplot(3, 4, i);
    try
        bode(system_links{i, 1}, omega);
        grid on;
        title(system_links{i, 2}, 'FontSize', 8);
    catch
        text(0.2, 0.5, 'АЧХ=const', 'FontSize', 8);
        grid on;
        title(system_links{i, 2}, 'FontSize', 8);
    end
end
figure('Name', 'АФЧХ всех звеньев', 'NumberTitle', 'off');
for i = 1:num_links
    subplot(3, 4, i);
    try
        nyquist(system_links{i, 1}, omega);
        grid on;
        title(system_links{i, 2}, 'FontSize', 8);
    catch
        text(0.2, 0.5, 'Точка (2.4, 0)', 'FontSize', 8);
        grid on;
        title(system_links{i, 2}, 'FontSize', 8);
    end
end

run_sensitivity_analysis('Влияние коэффициента K', [0.5, 1.0, 2.4, 5.0, 10.0], true, num_links, K, T, T1, T2, T1_5, s, omega, system_links);
run_sensitivity_analysis('Влияние постоянной времени T на ЛАЧХ всех звеньев', [0.1, 0.4, 1.0, 2.5, 5.0], false, num_links, K, T, T1, T2, T1_5, s, omega, system_links);

function W = get_link_formula(idx, curr_K, curr_T, curr_T1, curr_T2, curr_T1_5, s)
    switch idx
        case 1, W = curr_K / s;
        case 2, W = curr_K * s;
        case 3, W = tf(curr_K, 1);
        case 4, W = curr_K / (curr_T * s + 1);
        case 5, W = curr_K / (curr_T2^2 * s^2 + curr_T1_5 * s + 1);
        case 6, W = curr_K / (curr_T2^2 * s^2 + curr_T1 * s + 1);
        case 7, W = curr_K / (curr_T^2 * s^2 + 1);
        case 8, W = curr_K / (s * (curr_T * s + 1));
        case 9, W = (curr_K * s) / (curr_T * s + 1);
        case 10, W = curr_K * (curr_T * s + 1);
        case 11, W = (curr_K * (curr_T * s + 1)) / s;
    end
end

%Проверка влияния коэффициентов K и T на нашу систему
function run_sensitivity_analysis(fig_name, values, is_K_analysis, num_links, K, T, T1, T2, T1_5, s, omega, system_links)
    colors = lines(length(values));
    figure('Name', fig_name, 'NumberTitle', 'off');
    for i = 1:num_links
        subplot(3, 4, i);
        hold on;
        labels = cell(1, length(values)); % Пустая легенда на графике с перечислением параметров
        for idx = 1:length(values)
            val = values(idx);
            if is_K_analysis
                W_temp = get_link_formula(i, val, T, T1, T2, T1_5, s);
                label_str = sprintf('K=%.1f', val);
            else
                W_temp = get_link_formula(i, K, val, T1, T2, val, s);
                label_str = sprintf('T=%.1f', val);
            end
            [mag, ~, wout] = bode(W_temp, omega); % Вычисляет амплитуду и частоту для текущего звена
            mag = squeeze(mag);
            if isscalar(mag), mag = mag * ones(size(wout)); end 
            semilogx(wout, 20*log10(mag), 'Color', colors(idx, :), 'LineWidth', 1.1);
            labels{idx} = label_str;
        end
        grid on;
        title(system_links{i, 2}, 'FontSize', 8);
        xlabel('\omega, rad/s', 'FontSize', 7);
        ylabel('L(\omega), dB', 'FontSize', 7);
        if is_K_analysis && i == 1
            legend(labels, 'FontSize', 5, 'Location', 'best');
        elseif ~is_K_analysis && i == 1
            legend(labels, 'FontSize', 5, 'Location', 'best');
        end
    end
end
