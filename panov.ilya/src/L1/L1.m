syms s t;
K = 1.1;
T = 0.2;
T1 = 0.3;
T2 = 0.5;

W = {K/s, K*s, K, K/(T*s+1), K/(T2^2*s^2+T1*s+1), K/(T2^2*s^2+T1*s+1), K/(T*s^2+1), K/(s*(T*s+1)), K*s/(T*s+1), K*(T*s+1), K*(T*s+1)/s};

time = 0:0.01:10;
figure_h = figure('Name', 'Ступенчатое воздействие', 'NumberTitle', 'off', 'Position', [50 50 1200 800]);
figure_w = figure('Name', 'Импульсное воздействие', 'NumberTitle', 'off', 'Position', [80 80 1200 800]);

for i = 1:numel(W)
    U_step = laplace(heaviside(t));

    Y_step = W{i} * U_step;

    h = ilaplace(Y_step, s, t);

    U_impulse = laplace(dirac(t));

    Y_impulse = W{i} * U_impulse;

    w = ilaplace(Y_impulse, s, t);

    figure(figure_h);
    subplot(3, 4, i);
    hold on;
    grid on;

    h_plot = subs(h, dirac(1, t), 0);
    h_plot = subs(h_plot, dirac(t), 0);
    y = double(subs(h_plot, t, time)) + zeros(size(time));
    plot(time, y);

    h_delta = simplify(h - h_plot);
    if ~isequal(h_delta, sym(0))
        text(0.03, 0.85, ['$' latex(h_delta) '$'], ...
            'Units', 'normalized', 'Interpreter', 'latex');
    end

    title(['W_', num2str(i)]);
    xlabel('Время t, с');
    ylabel('h(t)');
    hold off;

    figure(figure_w);
    subplot(3, 4, i);
    hold on;
    grid on;

    w_plot = subs(w, dirac(1, t), 0);
    w_plot = subs(w_plot, dirac(t), 0);
    y = double(subs(w_plot, t, time)) + zeros(size(time));
    plot(time, y);

    w_delta = simplify(w - w_plot);
    if ~isequal(w_delta, sym(0))
        text(0.03, 0.85, ['$' latex(w_delta) '$'], ...
            'Units', 'normalized', 'Interpreter', 'latex');
    end

    title(['W_', num2str(i)]);
    xlabel('Время t, с');
    ylabel('w(t)');
    hold off;
end

K_values = 1:10;
figure_extra = figure('Name', 'Доп', 'NumberTitle', 'off', 'Position', [100 100 1200 800]);

subplot(2, 2, 1);
hold on;
grid on;

for k = K_values
    W_current = k/(T2^2*s^2 + T1*s + 1);
    h_k = ilaplace(W_current*U_step, s, t);
    y = double(subs(h_k, t, time));
    plot(time, y, 'DisplayName', ['K = ', num2str(k)]);
end

title('W_6: влияние K');
xlabel('Время t, с');
ylabel('h(t)');
legend('show', 'Location', 'northwest', 'NumColumns', 2);
hold off;

T_values = T*[1/5 1/2 1 2 5];

subplot(2, 2, 2);
hold on;
grid on;

for tau = T_values
    W_current = K/(tau*s + 1);
    h_k = ilaplace(W_current*U_step, s, t);
    y = double(subs(h_k, t, time));
    plot(time, y, 'DisplayName', ['T = ', num2str(tau)]);
end

title('W_4: влияние T');
xlabel('Время t, с');
ylabel('h(t)');
legend('show', 'Location', 'best');
hold off;

T1_values = T1*[1 1/2 1/5 1/10];

subplot(2, 2, 3);
hold on;
grid on;

for tau1 = T1_values
    W_current = K/(T2^2*s^2 + tau1*s + 1);
    h_k = ilaplace(W_current*U_step, s, t);
    y = double(subs(h_k, t, time));
    plot(time, y, 'DisplayName', ['T1 = ', num2str(tau1)]);
end

title('W_6: влияние T_1');
xlabel('Время t, с');
ylabel('h(t)');
legend('show', 'Location', 'best');
hold off;

s_tf = tf('s');
W_check = K/(T*s_tf + 1);
W_current = K/(T*s + 1);
h_check = ilaplace(W_current*U_step, s, t);
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
