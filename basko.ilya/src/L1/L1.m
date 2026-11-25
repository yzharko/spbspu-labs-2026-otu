clear; clc; close all;

K0 = 0.9;
T = 0.6;
T1 = 0.8;
T2 = 0.7;
T1a = 1.6; 

K_values = K0 * [0.5 1 2];
time = (0:0.01:10)';
colors = lines(numel(K_values));

link_names = {
    'Интегрирующее'
    'Идеальное дифференцирующее'
    'Усилительное'
    'Апериодическое 1-го порядка'
    'Апериодическое 2-го порядка'
    'Колебательное'
    'Консервативное'
    'Инерционное интегрирующее'
    'Реальное дифференцирующее'
    'Форсирующее 1-го порядка'
    'Изодромное'
};
link_count = numel(link_names);

syms s t K;
W = {
    K/s
    K*s
    K
    K/(T*s + 1)
    K/(T2^2*s^2 + T1a*s + 1)
    K/(T2^2*s^2 + T1*s + 1)
    K/(T*s^2 + 1)
    K/(s*(T*s + 1))
    K*s/(T*s + 1)
    K*(T*s + 1)
    K*(T*s + 1)/s
};

U_step = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);
h = cell(link_count, 1);
w = cell(link_count, 1);
h_regular = cell(link_count, 1);
w_regular = cell(link_count, 1);
h_singular = cell(link_count, 1);
w_singular = cell(link_count, 1);

fprintf('K = %.1f; T = %.1f; T1 = %.1f; T2 = %.1f\n', K0, T, T1, T2);
fprintf('Для апериодического звена 5 принято T1 = %.1f.\n', T1a);
fprintf('U_step(s) = %s; U_impulse(s) = %s\n', char(U_step), char(U_impulse));

for i = 1:link_count
    Y_step = simplify(W{i} * U_step);
    Y_impulse = simplify(W{i} * U_impulse);
    h{i} = simplify(ilaplace(Y_step, s, t));
    w{i} = simplify(ilaplace(Y_impulse, s, t));

    h_regular{i} = simplify(subs(h{i}, dirac(1,t), 0));
    h_regular{i} = simplify(subs(h_regular{i}, dirac(t), 0));
    w_regular{i} = simplify(subs(w{i}, dirac(1,t), 0));
    w_regular{i} = simplify(subs(w_regular{i}, dirac(t), 0));
    h_singular{i} = simplify(subs(h{i} - h_regular{i}, K, K0));
    w_singular{i} = simplify(subs(w{i} - w_regular{i}, K, K0));

    fprintf('\n%d. %s\n', i, link_names{i});
    fprintf('W(s) = %s\n', char(subs(W{i}, K, K0)));
    fprintf('Y_step(s) = %s\nh(t) = %s\n', ...
        char(subs(Y_step, K, K0)), char(subs(h{i}, K, K0)));
    fprintf('Y_impulse(s) = %s\nw(t) = %s\n', ...
        char(subs(Y_impulse, K, K0)), char(subs(w{i}, K, K0)));
end

s_tf = tf('s');
models = {
    K0/s_tf
    K0*s_tf
    tf(K0)
    K0/(T*s_tf + 1)
    K0/(T2^2*s_tf^2 + T1a*s_tf + 1)
    K0/(T2^2*s_tf^2 + T1*s_tf + 1)
    K0/(T*s_tf^2 + 1)
    K0/(s_tf*(T*s_tf + 1))
    K0*s_tf/(T*s_tf + 1)
    K0*(T*s_tf + 1)
    K0*(T*s_tf + 1)/s_tf
};
h_check = cell(link_count, 1);
w_check = cell(link_count, 1);

fprintf('\nПроверка обычных частей откликов при K = %.1f:\n', K0);
for i = 1:link_count
    if i == 2
        h_check{i} = zeros(size(time));
    elseif i == 10
        h_check{i} = K0 * ones(size(time));
    else
        h_check{i} = step(models{i}, time);
    end

    switch i
        case {2, 3, 10}
            w_check{i} = zeros(size(time));
        case 9
            w_check{i} = impulse(-K0/T/(T*s_tf + 1), time);
        case 11
            w_check{i} = impulse(K0/s_tf, time);
        otherwise
            w_check{i} = impulse(models{i}, time);
    end

    h0 = subs(h_regular{i}, K, K0);
    w0 = subs(w_regular{i}, K, K0);
    yh = double(subs(h0, t, time)) + zeros(size(time));
    yw = double(subs(w0, t, time)) + zeros(size(time));
    fprintf('%2d. max|h - h_check| = %.3g; max|w - w_check| = %.3g\n', ...
        i, max(abs(yh - h_check{i}(:))), max(abs(yw - w_check{i}(:))));

    if i == 2 || i == 10
        rh = simplify(laplace(h{i}, t, s) - W{i}*U_step);
        rw = simplify(laplace(w{i}, t, s) - W{i}*U_impulse);
        fprintf('    Проверка полных откликов по Лапласу: h = %s; w = %s\n', ...
            char(rh), char(rw));
    end
end

plot_functions = {h_regular, w_regular};
plot_checks = {h_check, w_check};
plot_singular = {h_singular, w_singular};
window_names = {'Переходные характеристики: обычные части', ...
                'Импульсные характеристики: обычные части'};
y_labels = {'h_{reg}(t)', 'w_{reg}(t)'};
legend_labels = [compose("K = %.2g", K_values), ...
                 compose("Проверка, K = %.2g", K0)];

for kind = 1:2
    current_formulas = plot_functions{kind};
    checks = plot_checks{kind};
    singular = plot_singular{kind};

    figure('Name', window_names{kind}, 'NumberTitle', 'off', ...
        'Color', 'w', 'Position', [60 60 1280 800]);
    layout = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    curves = gobjects(1, numel(K_values) + 1);

    for i = 1:link_count
        ax = nexttile(layout);
        hold(ax, 'on'); grid(ax, 'on');

        for j = 1:numel(K_values)
            current_function = subs(current_formulas{i}, K, K_values(j));
            y = double(subs(current_function, t, time)) + zeros(size(time));
            curves(j) = plot(ax, time, y, 'Color', colors(j,:), 'LineWidth', 1.3);
        end

        points = 1:40:numel(time);
        curves(end) = plot(ax, time(points), checks{i}(points), ...
            'ko', 'MarkerSize', 3);
        if ~isequal(singular{i}, sym(0))
            text(ax, 0.03, 0.94, ['$K = ' num2str(K0) ':\quad ' latex(singular{i}) '$'], ...
                'Units', 'normalized', 'Interpreter', 'latex', ...
                'VerticalAlignment', 'top', 'BackgroundColor', 'w');
        end

        title(ax, sprintf('%d. %s', i, link_names{i}), ...
            'Interpreter', 'none', 'FontSize', 9);
        xlabel(ax, 't, с'); ylabel(ax, y_labels{kind});
        xlim(ax, [time(1) time(end)]);
    end

    lgd = legend(curves, legend_labels, 'Location', 'northwest');
    lgd.Layout.Tile = 12;
    title(layout, window_names{kind});
end

T_values = T * [1/5 1/2 1 2 5];
time_T = (0:0.01:16)';
figure('Name', 'Влияние T', 'NumberTitle', 'off', 'Color', 'w');
hold on; grid on;

fprintf('\nИзменение T: время достижения 95%% конечного значения\n');
fprintf('T, с    h''(0+)    t95, с\n');
for current_T = T_values
    model = K0/(current_T*s_tf + 1);
    [y, tout] = step(model, time_T);
    plot(tout, y, 'LineWidth', 1.4, ...
        'DisplayName', sprintf('T = %.2g с', current_T));
    fprintf('%.2f    %.3f     %.3f\n', ...
        current_T, K0/current_T, -current_T*log(0.05));
end

title('Апериодическое звено: влияние T');
xlabel('t, с'); ylabel('h(t)');
legend('Location', 'southeast');

T1_values = T1 * [1 1/2 1/5 1/10 0];
time_damping = (0:0.01:40)';
figure('Name', 'Влияние T1', 'NumberTitle', 'off', 'Color', 'w');
hold on; grid on;

fprintf('\nT1, с    Перерегулирование, %%    Постоянная затухания, с    Оценка t5%%, с\n');
for current_T1 = T1_values
    model = K0/(T2^2*s_tf^2 + current_T1*s_tf + 1);
    [y, tout] = step(model, time_damping);
    plot(tout, y, 'LineWidth', 1.3, ...
        'DisplayName', sprintf('T1 = %.2g с', current_T1));

    alpha = current_T1/(2*T2^2);
    omega = sqrt(4*T2^2 - current_T1^2)/(2*T2^2);
    overshoot = 100 * exp(-pi*alpha/omega);
    if current_T1 > 0
        decay_time = 1/alpha;
        settling_bound = log(20/sqrt(1 - (current_T1/(2*T2))^2))/alpha;
    else
        decay_time = Inf;
        settling_bound = Inf;
    end
    fprintf('%.2f       %.2f                    %.3f                    %.3f\n', ...
        current_T1, overshoot, decay_time, settling_bound);
end

title(sprintf('Колебательное звено: влияние T1 при T2 = %.2g с', T2));
xlabel('t, с'); ylabel('h(t)');
legend('Location', 'northeast');
