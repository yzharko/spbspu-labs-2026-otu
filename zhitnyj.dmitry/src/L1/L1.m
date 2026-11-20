clear; clc; close all;

K = 1.5;
T = 0.4;
T1 = 0.6;
T2 = 0.8;
T1_5 = 2.0; % Для звена №5: T1_5 > 2*T2.
syms s t;

names = {
    'Интегрирующее';
    'Идеальное дифференцирующее';
    'Усилительное';
    'Апериодическое 1-го порядка';
    'Апериодическое 2-го порядка';
    'Колебательное';
    'Консервативное';
    'Инерционное интегрирующее';
    'Реальное дифференцирующее';
    'Форсирующее 1-го порядка';
    'Изодромное'
};

W = {
    K/s;
    K*s;
    sym(K);
    K/(T*s + 1);
    K/(sym(T2, 'r')^2*s^2 + T1_5*s + 1);
    K/(sym(T2, 'r')^2*s^2 + T1*s + 1);
    K/(T*s^2 + 1);
    K/(s*(T*s + 1));
    K*s/(T*s + 1);
    K*(T*s + 1);
    K*(T*s + 1)/s
};

U_step = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);
time = 0:0.01:15;
h = cell(size(W));
w = cell(size(W));
h_analytic = zeros(numel(time), numel(W));
w_analytic = zeros(numel(time), numel(W));
h_axes = gobjects(size(W));
w_axes = gobjects(size(W));

figure_h = figure('Name', 'Переходные характеристики', ...
    'NumberTitle', 'off', 'Theme', 'light', 'Color', 'white', 'Position', [50 50 1400 850]);
figure_w = figure('Name', 'Импульсные характеристики', ...
    'NumberTitle', 'off', 'Theme', 'light', 'Color', 'white', 'Position', [80 80 1400 850]);

for i = 1:numel(W)
    Y_step = W{i}*U_step;
    Y_impulse = W{i}*U_impulse;
    h{i} = simplify(ilaplace(Y_step, s, t));
    w{i} = simplify(ilaplace(Y_impulse, s, t));

    h_regular = subs(h{i}, dirac(1, t), 0);
    h_regular = subs(h_regular, dirac(t), 0);
    h_delta = simplify(h{i} - h_regular);
    h_regular = subs(h_regular, heaviside(t), 1);
    h_values = double(subs(h_regular, t, time)) + zeros(size(time));
    h_values(1) = double(limit(h_regular, t, 0, 'right'));
    h_analytic(:,i) = h_values(:);

    w_regular = subs(w{i}, dirac(1, t), 0);
    w_regular = subs(w_regular, dirac(t), 0);
    w_delta = simplify(w{i} - w_regular);
    w_regular = subs(w_regular, heaviside(t), 1);
    w_values = double(subs(w_regular, t, time)) + zeros(size(time));
    w_values(1) = double(limit(w_regular, t, 0, 'right'));
    w_analytic(:,i) = w_values(:);

    fprintf('\n%d. %s\nW(s) = %s\nh(t) = %s\nw(t) = %s\n', ...
        i, names{i}, char(W{i}), char(h{i}), char(w{i}));
    figure(figure_h);
    h_axes(i) = subplot(3, 4, i);
    plot(time, h_values, 'b-');
    grid on;
    title(sprintf('%d. %s', i, names{i}));
    xlabel('t, с'); ylabel('h(t)');
    if ~isequal(h_delta, sym(0))
        text(0.03, 0.85, ['$+ ' latex(h_delta) '$'], ...
            'Units', 'normalized', 'Interpreter', 'latex');
    end

    figure(figure_w);
    w_axes(i) = subplot(3, 4, i);
    plot(time, w_values, 'b-');
    grid on;
    title(sprintf('%d. %s', i, names{i}));
    xlabel('t, с'); ylabel('w(t)');
    if ~isequal(w_delta, sym(0))
        text(0.03, 0.85, ['$+ ' latex(w_delta) '$'], ...
            'Units', 'normalized', 'Interpreter', 'latex');
    end
end

fprintf('\nПроверка через Control System Toolbox\n');
G = {
    tf(K, [1 0]);
    tf([K 0], 1);
    tf(K, 1);
    tf(K, [T 1]);
    tf(K, [T2^2 T1_5 1]);
    tf(K, [T2^2 T1 1]);
    tf(K, [T 0 1]);
    tf(K, [T 1 0]);
    tf([K 0], [T 1]);
    tf([K*T K], 1);
    tf([K*T K], [1 0])
};
check_tolerance = 1e-7;
step_error = NaN(size(W));
impulse_error = NaN(size(W));
fprintf('Допуск: %.0e\n', check_tolerance);
for i = 1:numel(W)
    if ~isproper(G{i})
        fprintf('Звено №%d: step и impulse неприменимы.\n', i);
        continue;
    end
    [h_check, time_h] = step(G{i}, time);
    [w_check, time_w] = impulse(G{i}, time);
    step_error(i) = max(abs(h_analytic(:,i) - h_check(:)));
    impulse_error(i) = max(abs(w_analytic(:,i) - w_check(:)));
    fprintf('Звено №%d: ошибка step %.3g; ошибка impulse %.3g\n', ...
        i, step_error(i), impulse_error(i));
    assert(step_error(i) < check_tolerance && impulse_error(i) < check_tolerance, ...
        ['Проверка не пройдена для звена ', num2str(i)]);
    hold(h_axes(i), 'on');
    plot(h_axes(i), time_h, h_check, 'r--');
    hold(h_axes(i), 'off');
    hold(w_axes(i), 'on');
    plot(w_axes(i), time_w, w_check, 'r--');
    hold(w_axes(i), 'off');
    if i == 1
        legend(h_axes(i), 'Аналитически', 'step', 'Location', 'best');
        legend(w_axes(i), 'Аналитически', 'impulse', 'Location', 'best');
    end
end

K_values = K*[0.5 1 2];
figure('Name', 'Влияние коэффициента усиления на временные характеристики', ...
    'NumberTitle', 'off', 'Theme', 'light', 'Color', 'white');
for k = K_values
    W4 = tf(k, [T 1]);
    W6 = tf(k, [T2^2 T1 1]);
    systems = {W4, W6};
    for j = 1:2
        subplot(2, 2, j);
        hold on; grid on;
        plot(time, step(systems{j}, time), ...
            'DisplayName', sprintf('K = %g', k));
        title(names{2*j + 2});
        set(gca, 'UserData', 2*j + 2);
        xlabel('t, с'); ylabel('h(t)');
        subplot(2, 2, j + 2);
        hold on; grid on;
        plot(time, impulse(systems{j}, time), ...
            'DisplayName', sprintf('K = %g', k));
        title(names{2*j + 2});
        set(gca, 'UserData', 2*j + 2);
        xlabel('t, с'); ylabel('w(t)');
    end
end
for j = 1:4
    subplot(2, 2, j);
    legend('show', 'Location', 'best');
    hold off;
end

T_values = T*[1/5 1/2 1 2 5];
time_T = 0:0.01:30;
figure('Name', 'Влияние постоянной времени на характеристики апериодического звена', ...
    'NumberTitle', 'off', 'Theme', 'light', 'Color', 'white');
fprintf('\nЗвено №4, K = %g; время установления по полосе 5%%\n', K);
for tau = T_values
    W4 = tf(K, [tau 1]);
    subplot(2, 1, 1);
    hold on; grid on;
    plot(time_T, step(W4, time_T), ...
        'DisplayName', sprintf('T = %g', tau));
    title('Переходная характеристика при различных T');
    xlabel('t, с'); ylabel('h(t)');
    subplot(2, 1, 2);
    hold on; grid on;
    plot(time_T, impulse(W4, time_T), ...
        'DisplayName', sprintf('T = %g', tau));
    title('Импульсная характеристика при различных T');
    xlabel('t, с'); ylabel('w(t)');
    fprintf('T = %.2f с; t_уст = %.4f с; w(0+) = %.4f\n', ...
        tau, -tau*log(0.05), K/tau);
end
for j = 1:2
    subplot(2, 1, j);
    legend('show', 'Location', 'best');
    hold off;
end

T1_values = [T1 T1/2 T1/5 T1/10 0];
time_T1 = 0:0.01:160;
figure('Name', 'Влияние коэффициента T1 на характеристики колебательного звена', ...
    'NumberTitle', 'off', 'Theme', 'light', 'Color', 'white');
fprintf('\nЗвено №6, K = %g, T2 = %g; полоса установления 5%%\n', K, T2);
for tau1 = T1_values
    W6 = tf(K, [T2^2 tau1 1]);
    y_h = step(W6, time_T1);
    y_w = impulse(W6, time_T1);
    alpha = tau1/(2*T2^2);
    zeta = tau1/(2*T2);
    omega = sqrt(1 - zeta^2)/T2;
    overshoot = 100*exp(-pi*zeta/sqrt(1 - zeta^2));
    last_outside = find(abs(y_h - K) > 0.05*abs(K), 1, 'last');
    if tau1 == 0 || (~isempty(last_outside) && last_outside == numel(time_T1))
        settling_time = Inf;
    elseif isempty(last_outside)
        settling_time = 0;
    else
        settling_time = time_T1(last_outside + 1);
    end
    fprintf(['T1 = %.2f; zeta = %.4f; alpha = %.6f; omega = %.6f; ' ...
        'перерегулирование = %.2f%%; t_уст = %.2f с\n'], ...
        tau1, zeta, alpha, omega, overshoot, settling_time);
    subplot(2, 1, 1);
    hold on; grid on;
    plot(time_T1, y_h, 'DisplayName', sprintf('T1 = %g', tau1));
    title('Переходная характеристика при различных T_1');
    xlabel('t, с'); ylabel('h(t)');
    subplot(2, 1, 2);
    hold on; grid on;
    plot(time_T1, y_w, 'DisplayName', sprintf('T1 = %g', tau1));
    title('Импульсная характеристика при различных T_1');
    xlabel('t, с'); ylabel('w(t)');
end
for j = 1:2
    subplot(2, 1, j);
    legend('show', 'Location', 'best');
    hold off;
end
