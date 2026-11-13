clear; close all; clc;
set(groot, 'defaultAxesColor', 'w', 'defaultAxesXColor', 'k', ...
    'defaultAxesYColor', 'k', 'defaultTextColor', 'k', ...
    'defaultAxesGridColor', [0.75 0.79 0.83]);

K = 2.5;
T = 0.8;
T1 = 0.8;
T2 = 1.0;

output_dir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
t = linspace(0, 12, 2401);

names = { ...
    'Идеальное интегрирующее', 'Идеальное дифференцирующее', ...
    'Усилительное', 'Апериодическое 1-го порядка', ...
    'Апериодическое 2-го порядка', 'Колебательное', ...
    'Консервативное', 'Инерционное интегрирующее', ...
    'Реальное дифференцирующее', 'Форсирующее 1-го порядка', ...
    'Изодромное'};

h = cell(1, 11);
w = cell(1, 11);
step_singular = cell(1, 11);
impulse_singular = cell(1, 11);

h{1} = K*t;
w{1} = K + zeros(size(t));

h{2} = zeros(size(t));
w{2} = zeros(size(t));
step_singular{2} = 'K delta(t)';
impulse_singular{2} = 'K delta''(t)';

h{3} = K + zeros(size(t));
w{3} = zeros(size(t));
impulse_singular{3} = 'K delta(t)';

h{4} = K*(1-exp(-t/T));
w{4} = (K/T)*exp(-t/T);

if T1 >= 2*T2
    [h{5}, w{5}] = second_order_response(t, K, T1, T2);
else
    h{5} = nan(size(t));
    w{5} = nan(size(t));
end

[h{6}, w{6}] = second_order_response(t, K, T1, T2);

omega0 = 1/sqrt(T);
h{7} = K*(1-cos(omega0*t));
w{7} = K*omega0*sin(omega0*t);

h{8} = K*(t-T*(1-exp(-t/T)));
w{8} = K*(1-exp(-t/T));

h{9} = (K/T)*exp(-t/T);
w{9} = -(K/T^2)*exp(-t/T);
impulse_singular{9} = '(K/T) delta(t)';

h{10} = K + zeros(size(t));
w{10} = zeros(size(t));
step_singular{10} = 'K*T delta(t)';
impulse_singular{10} = 'K delta(t) + K*T delta''(t)';

h{11} = K*T + K*t;
w{11} = K + zeros(size(t));
impulse_singular{11} = 'K*T delta(t)';

plot_collection(t, h, names, step_singular, ...
    'Переходные характеристики h(t)', 'h(t)', ...
    fullfile(output_dir, 'step_responses.png'));
plot_collection(t, w, names, impulse_singular, ...
    'Импульсные характеристики w(t): регулярная часть и сингулярные слагаемые', ...
    'w(t)', fullfile(output_dir, 'impulse_responses.png'));

K_values = [K/2 K 2*K];
f = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 720]);
tl = tiledlayout(2, 1, 'Padding', 'compact');
title(tl, 'Влияние K: исходный вариант и масштабирование', 'Color', 'k');
nexttile; hold on; grid on;
for j = 1:numel(K_values)
    Kj = K_values(j);
    plot(t, Kj*(1-exp(-t/T)), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('K = %.2f', Kj));
end
title('Апериодическое звено №4', 'Color', 'k');
xlabel('t, с'); ylabel('h(t)');
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
legend('Location', 'southeast', 'Color', 'w', 'TextColor', 'k');
nexttile; hold on; grid on;
for j = 1:numel(K_values)
    Kj = K_values(j);
    [hj, ~] = second_order_response(t, Kj, T1, T2);
    plot(t, hj, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('K = %.2f', Kj));
end
title('Колебательное звено №6', 'Color', 'k');
xlabel('t, с'); ylabel('h(t)');
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
legend('Location', 'northeast', 'Color', 'w', 'TextColor', 'k');
exportgraphics(f, fullfile(output_dir, 'gain_variation.png'), 'Resolution', 180);
close(f);

T_values = [T/5, T/2, T, 2*T, 5*T];
t_slow = linspace(0, 20, 3001);
f = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 520]);
hold on; grid on;
for j = 1:numel(T_values)
    Tj = T_values(j);
    plot(t_slow, K*(1-exp(-t_slow/Tj)), 'LineWidth', 1.6, ...
        'DisplayName', sprintf('T = %.2f с; t_{95} = %.2f с', Tj, -Tj*log(0.05)));
end
xlabel('Время t, с'); ylabel('h(t)');
title('Звено №4: влияние постоянной времени', 'Color', 'k');
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
legend('Location', 'southeast', 'Color', 'w', 'TextColor', 'k');
exportgraphics(f, fullfile(output_dir, 'time_variation.png'), 'Resolution', 180);
close(f);

T1_values = [0.8 0.4 0.2];
t_damp = linspace(0, 55, 5501);
f = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 520]);
hold on; grid on;
fprintf('Звено №6, T2 = %.2f с:\n', T2);
fprintf(' T1, с | zeta | перерегулирование, %% | оценка t_2%%, с\n');
for j = 1:numel(T1_values)
    T1j = T1_values(j);
    [hj, ~] = second_order_response(t_damp, K, T1j, T2);
    zeta = T1j/(2*T2);
    overshoot = 100*exp(-pi*zeta/sqrt(1-zeta^2));
    settling_approx = -log(0.02)*2*T2^2/T1j;
    fprintf(' %5.2f | %.2f | %19.2f | %16.2f\n', ...
        T1j, zeta, overshoot, settling_approx);
    plot(t_damp, hj, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('T_1=%.1f с; M_p=%.1f%%', T1j, overshoot));
end
yline(K, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
xlabel('Время t, с'); ylabel('h(t)');
title('Звено №6: уменьшение T_1 при неизменном T_2', 'Color', 'k');
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
legend('Location', 'northeast', 'Color', 'w', 'TextColor', 'k');
exportgraphics(f, fullfile(output_dir, 'damping_variation.png'), 'Resolution', 180);
close(f);

fprintf('Графики сохранены в %s\n', output_dir);

function [h, w] = second_order_response(t, K, T1, T2)
    roots_s = roots([T2^2 T1 1]);
    r1 = roots_s(1);
    r2 = roots_s(2);
    if abs(r1-r2) < 1e-10
        h = K*(1-exp(r1*t).*(1-r1*t));
        w = (K/T2^2)*t.*exp(r1*t);
    else
        h = K*(1+(r2*exp(r1*t)-r1*exp(r2*t))/(r1-r2));
        w = (K/T2^2)*(exp(r1*t)-exp(r2*t))/(r1-r2);
    end
    h = real(h);
    w = real(w);
end

function plot_collection(t, curves, names, singular, heading, y_label, file_name)
    f = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1500 1100]);
    tl = tiledlayout(4, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, heading, 'FontSize', 17, 'FontWeight', 'bold', 'Color', 'k');
    for i = 1:numel(curves)
        nexttile;
        if all(isnan(curves{i}))
            text(0.5, 0.5, 'Условие T_1 >= 2T_2 не выполнено; график не строится', ...
                'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                'FontSize', 9, 'Color', 'k');
            axis off;
            title(sprintf('%d. %s', i, names{i}), 'FontSize', 10, 'Color', 'k');
            continue;
        end
        plot(t, curves{i}, 'b', 'LineWidth', 1.5);
        grid on; xlim([0 t(end)]);
        set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
        title(sprintf('%d. %s', i, names{i}), 'FontSize', 10, 'Color', 'k');
        xlabel('t, с'); ylabel(y_label);
        if max(abs(curves{i})) < 1e-10
            ylim([-1 1]);
        end
        if ~isempty(singular{i})
            hold on;
            xline(0, 'k--', 'LineWidth', 1.4);
            yl = ylim;
            text(t(end)*0.04, yl(2)-0.12*(yl(2)-yl(1)), ...
                ['Слагаемое: ' singular{i}], 'Color', 'k', ...
                'FontSize', 9, 'Interpreter', 'none');
        end
    end
    exportgraphics(f, file_name, 'Resolution', 180);
    close(f);
end
