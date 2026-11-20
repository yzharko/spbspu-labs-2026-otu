close all; clear; clc;

pkg load symbolic
pkg load control

syms t s

K = 2.6;
T = 0.3;
T1 = 0.6;
T2 = 0.5;
T1_5 = 1.2;

W = {
    K / s;
    K * s;
    K;
    K / (T*s + 1);
    K / (T2^2*s^2 + T1_5*s + 1);
    K / (T2^2*s^2 + T1*s + 1);
    K / (T*s^2 + 1);
    K / (s*(T*s + 1));
    K*s / (T*s + 1);
    K*(T*s + 1);
    K*(T*s + 1) / s
};

u_step = heaviside(t);
U_step = laplace(u_step, t, s);

u_impulse = dirac(t);
U_impulse = laplace(u_impulse, t, s);

Y_step = cell(size(W));
Y_impulse = cell(size(W));

for i = 1:numel(W)
    Y_step{i} = simplify(W{i} * U_step);
    Y_impulse{i} = simplify(W{i} * U_impulse);
end

h = cell(size(W));
w = cell(size(W));

for i = 1:numel(W)
    h{i} = simplify(ilaplace(Y_step{i}, s, t));
    w{i} = simplify(ilaplace(Y_impulse{i}, s, t));
end

time = 0.01:0.01:5;
delta = dirac(t);
delta_derivative = diff(delta, t);

for i = 1:numel(W)
    h_expr = h{i};
    w_expr = w{i};

    h_expr = subs(h_expr, delta_derivative, sym(0));
    w_expr = subs(w_expr, delta_derivative, sym(0));

    h_expr = subs(h_expr, delta, sym(0));
    w_expr = subs(w_expr, delta, sym(0));

    h_expr = subs(h_expr, heaviside(t), sym(1));
    w_expr = subs(w_expr, heaviside(t), sym(1));

    h_func = function_handle(h_expr, 'vars', t);
    w_func = function_handle(w_expr, 'vars', t);

    h_values = h_func(time);
    w_values = w_func(time);

    if isscalar(h_values)
        h_values = h_values * ones(size(time));
    end

    if isscalar(w_values)
        w_values = w_values * ones(size(time));
    end

    figure(i);
    clf;

    subplot(2, 1, 1);
    plot(time, h_values, 'b', 'LineWidth', 2);
    grid on;
    title(sprintf('Звено %d: h(t)', i));
    xlabel('Время, с');
    ylabel('h(t)');

    subplot(2, 1, 2);
    plot(time, w_values, 'r', 'LineWidth', 2);
    grid on;
    title(sprintf('Звено %d: w(t)', i));
    xlabel('Время, с');
    ylabel('w(t)');
end

T_values = T * [1/5 1/2 1 2 5];
study_time = 0:0.01:10;

figure(12);
clf;
hold on;
grid on;

for Ti = T_values
    W_control = tf(K, [Ti 1]);
    [y, t_out] = step(W_control, study_time);

    plot(t_out, y, 'LineWidth', 1.5);
end

title('Влияние T на апериодическое звено №4');
xlabel('Время, с');
ylabel('h(t)');

labels = arrayfun(@(x) sprintf('T = %.2f', x), ...
                  T_values, 'UniformOutput', false);
legend(labels, 'Location', 'best');study_time = 0:0.01:20;
K_values = 1:10;

figure(13);
clf;
hold on;
grid on;

labels = cell(size(K_values));

for i = 1:numel(K_values)

    Ki = K_values(i);

    W_control = tf(Ki, [T2^2 T1 1]);

    [y, t_out] = step(W_control, study_time);

    plot(t_out, y, 'LineWidth', 1.5);

    labels{i} = sprintf('K = %d', Ki);
end

title('Влияние K на колебательное звено');
xlabel('Время, с');
ylabel('h(t)');
legend(labels, 'Location', 'eastoutside');

hold off;
hold off;

T1_values = T1 * [1 1/2 1/5 1/10];
study_time = 0:0.01:45;

figure(14);
clf;
hold on;
grid on;

for T1i = T1_values
    W_control = tf(K, [T2^2 T1i 1]);
    [y, t_out] = step(W_control, study_time);

    plot(t_out, y, 'LineWidth', 1.5);

    overshoot = (max(y) - K) / K * 100;

    idx = find(abs(y - K) > 0.02*K, 1, 'last');

    if isempty(idx)
        settle_time = 0;
    elseif idx == numel(t_out)
        settle_time = NaN;
    else
        settle_time = t_out(idx + 1);
    end

    fprintf('T1 = %.3f: перерегулирование = %.2f%%, ', ...
            T1i, overshoot);
    fprintf('время установления = %.2f с\n', settle_time);
end

title('Влияние T1 на колебательное звено №6');
xlabel('Время, с');
ylabel('h(t)');

labels = arrayfun(@(x) sprintf('T1 = %.3f', x), ...
                  T1_values, 'UniformOutput', false);
legend(labels, 'Location', 'best');
hold off;

check_time = [0 time];

for i = 1:numel(W)

    [num, den] = numden(sym(W{i}));

    W_control = tf(double(sym2poly(num)), ...
                   double(sym2poly(den)));

    fprintf('\nЗвено %d\n', i);

    if ~ismember(i, [2 10])
        h_expr = subs(h{i}, heaviside(t), sym(1));
        h_func = function_handle(h_expr, 'vars', t);
        h_analytic = h_func(time) + zeros(size(time));

        h_control = step(W_control, check_time);
        h_control = h_control(2:end);

        error_h = max(abs(h_analytic(:) - h_control(:)));

        fprintf('Ошибка h(t): %.6e\n', error_h);
    else
        disp('h(t): проверка пропущена');
    end

    if ~ismember(i, [2 3 9 10 11])

        w_expr = subs(w{i}, heaviside(t), sym(1));
        w_func = function_handle(w_expr, 'vars', t);
        w_analytic = w_func(time) + zeros(size(time));

        w_control = impulse(W_control, check_time);
        w_control = w_control(2:end);

        error_w = max(abs(w_analytic(:) - w_control(:)));

        fprintf('Ошибка w(t): %.6e\n', error_w);
    else
        disp('w(t): проверка пропущена');
    end

end

