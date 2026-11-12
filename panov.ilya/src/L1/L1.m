syms s t;
K = 1.1;
T = 0.2;
T1 = 0.3;
T2 = 0.5;
T1_5 = 1.2;

W = {
    K / s;
    K * s;
    K;
    K / (T * s + 1);
    K / (T2^2 * s^2 + T1_5 * s + 1);
    K / (T2^2 * s^2 + T1 * s + 1);
    K / (T * s^2 + 1);
    K / (s * (T * s + 1));
    K * s / (T * s + 1);
    K * (T * s + 1);
    K * (T * s + 1) / s
};

U_step = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);
time = 0:0.01:10;
step_error = zeros(size(W));
impulse_error = zeros(size(W));
delta_match = false(size(W));

figure_h = figure('Name', 'Ступенчатое воздействие', ...
    'NumberTitle', 'off', 'Position', [50 50 1200 800]);
figure_w = figure('Name', 'Импульсное воздействие', ...
    'NumberTitle', 'off', 'Position', [80 80 1200 800]);

for i = 1:numel(W)
    Y_step = W{i} * U_step;
    h = ilaplace(Y_step, s, t);

    Y_impulse = W{i} * U_impulse;
    w = ilaplace(Y_impulse, s, t);

    h_plot = subs(h, dirac(1, t), 0);
    h_plot = subs(h_plot, dirac(t), 0);
    h_delta = simplify(h - h_plot);
    h_values = double(subs(h_plot, t, time)) + zeros(size(time));
    h_values(1) = double(limit(h_plot, t, 0, 'right'));

    w_plot = subs(w, dirac(1, t), 0);
    w_plot = subs(w_plot, dirac(t), 0);
    w_delta = simplify(w - w_plot);
    w_values = double(subs(w_plot, t, time)) + zeros(size(time));
    w_values(1) = double(limit(w_plot, t, 0, 'right'));

    [num, den] = numden(sym(W{i}));
    num = sym2poly(num);
    den = sym2poly(den);
    [q, remainder] = deconv(num, den);
    q_0 = q(end);
    q_s = 0;
    if numel(q) == 2
        q_s = q(1);
    end

    W_regular = tf(remainder, den);
    [h_check, time_step] = step(W_regular + q_0, time);
    [w_check, time_impulse] = impulse(W_regular, time);
    h_delta_check = sym(q_s, 'r') * dirac(t);
    w_delta_check = sym(q_0, 'r') * dirac(t) ...
        + sym(q_s, 'r') * dirac(1, t);

    step_error(i) = max(abs(h_values(:) - h_check(:)));
    impulse_error(i) = max(abs(w_values(:) - w_check(:)));
    delta_match(i) = isAlways(simplify(h_delta - h_delta_check) == 0) ...
        && isAlways(simplify(w_delta - w_delta_check) == 0);
    assert(step_error(i) < 1e-8 && impulse_error(i) < 1e-8 ...
        && delta_match(i), ['Не совпали результаты для W_', num2str(i)]);

    figure(figure_h);
    subplot(3, 4, i);
    hold on;
    grid on;

    plot(time, h_values, 'b-', 'DisplayName', 'Аналитически');
    plot(time_step, h_check, 'r--', 'DisplayName', 'step');
    if ~isequal(h_delta, sym(0))
        text(0.03, 0.85, ['$' latex(h_delta) '$'], ...
            'Units', 'normalized', 'Interpreter', 'latex');
    end

    title(['W_', num2str(i)]);
    xlabel('Время t, с');
    ylabel('h(t)');
    if i == 1
        legend('show', 'Location', 'best');
    end
    hold off;

    figure(figure_w);
    subplot(3, 4, i);
    hold on;
    grid on;

    plot(time, w_values, 'b-', 'DisplayName', 'Аналитически');
    plot(time_impulse, w_check, 'r--', 'DisplayName', 'impulse');
    if ~isequal(w_delta, sym(0))
        text(0.03, 0.85, ['$' latex(w_delta) '$'], ...
            'Units', 'normalized', 'Interpreter', 'latex');
    end

    title(['W_', num2str(i)]);
    xlabel('Время t, с');
    ylabel('w(t)');
    if i == 1
        legend('show', 'Location', 'best');
    end
    hold off;
end

K_values = 1:10;
figure_extra = figure('Name', 'Доп', 'NumberTitle', 'off', ...
    'Position', [100 100 1200 800]);

subplot(2, 2, 1);
hold on;
grid on;

for k = K_values
    W_current = k / (T2^2 * s^2 + T1 * s + 1);
    h_k = ilaplace(W_current * U_step, s, t);
    y = double(subs(h_k, t, time));
    plot(time, y, 'DisplayName', ['K = ', num2str(k)]);
end

title('W_6: влияние K');
xlabel('Время t, с');
ylabel('h(t)');
legend('show', 'Location', 'northwest', 'NumColumns', 2);
hold off;

T_values = T * [1/5 1/2 1 2 5];

subplot(2, 2, 2);
hold on;
grid on;

for tau = T_values
    W_current = K / (tau * s + 1);
    h_k = ilaplace(W_current * U_step, s, t);
    y = double(subs(h_k, t, time));
    plot(time, y, 'DisplayName', ['T = ', num2str(tau)]);
end

title('W_4: влияние T');
xlabel('Время t, с');
ylabel('h(t)');
legend('show', 'Location', 'best');
hold off;

T1_values = T1 * [1 1/2 1/5 1/10];

subplot(2, 2, 3);
hold on;
grid on;

for tau1 = T1_values
    W_current = K / (T2^2 * s^2 + tau1 * s + 1);
    h_k = ilaplace(W_current * U_step, s, t);
    y = double(subs(h_k, t, time));
    plot(time, y, 'DisplayName', ['T1 = ', num2str(tau1)]);
end

title('W_6: влияние T_1');
xlabel('Время t, с');
ylabel('h(t)');
legend('show', 'Location', 'best');
hold off;

s_tf = tf('s');
W_check = K / (T * s_tf + 1);
W_current = K / (T * s + 1);
h_check = ilaplace(W_current * U_step, s, t);
y = double(subs(h_check, t, time));
[y_check, time_check] = step(W_check, time);

subplot(2, 2, 4);
hold on;
grid on;

plot(time, y, 'b-', 'DisplayName', 'Аналитически');
plot(time_check, y_check, 'r--', 'DisplayName', 'step');

title('W_4: проверка');
xlabel('Время t, с');
ylabel('h(t)');
legend('show', 'Location', 'best');
hold off;
