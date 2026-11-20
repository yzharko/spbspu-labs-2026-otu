syms s t K
K0 = 1.5;
T = 0.4;
T1 = 0.6;
T2 = 0.8;

K_values = K0*[1 2 3 4];
time = 0:0.01:10;
s_tf = tf('s');

names = {'Интегрирующее', ...
    'Идеальное дифференцирующее', ...
    'Усилительное', ...
    'Апериодическое 1-го порядка', ...
    'Апериодическое 2-го порядка', ...
    'Колебательное', ...
    'Консервативное', ...
    'Инерционное интегрирующее', ...
    'Реальное дифференцирующее', ...
    'Форсирующее', ...
    'Изодромное'};

W_all = [K/s; ...
    K*s; K; ...
    K/(T*s+1); ...
    K/(T2^2*s^2+T1*s+1); ...
    K/(T2^2*s^2+T1*s+1); ...
    K/(T*s^2+1); ...
    K/(s*(T*s+1)); ...
    K*s/(T*s+1); ...
    K*(T*s+1); ...
    K*(T*s+1)/s];

U_step = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);

h_all = cell(1, 11);
w_all = cell(1, 11);

figure(1); clf;
figure(2); clf;

for n = 1:11
    W = W_all(n);
    Y_step = W*U_step;
    h = simplify(ilaplace(Y_step, s, t));
    Y_impulse = W*U_impulse;
    w = simplify(ilaplace(Y_impulse, s, t));

    h_plot = subs(subs(h, dirac(1,t), 0), dirac(t), 0);
    w_plot = subs(subs(w, dirac(1,t), 0), dirac(t), 0);
    h_all{n} = h_plot;
    w_all{n} = w_plot;

    for p = 1:2
        figure(p);
        subplot(3, 4, n);
        hold on; grid on;

        for k = K_values
            if p == 1
                f = h_plot;
            else
                f = w_plot;
            end
            f_k = subs(f, K, k);
            y = double(subs(f_k, t, time)) + zeros(size(time));
            plot(time, y, 'DisplayName', ['K = ', num2str(k)]);
        end

        title(sprintf('%d. %s', n, names{n}), 'FontSize', 9);
        xlabel('Время t, с');
        if p == 1
            ylabel('h(t)');
        else
            ylabel('w(t)');
        end
        hold off;
    end
end

for p = 1:2
    figure(p);
    subplot(3, 4, 12); hold on;
    for k = K_values
        plot(NaN, NaN, 'DisplayName', ['K = ', num2str(k)]);
    end
    axis off;
    legend('show', 'Location', 'northwest');
end
figure(1); sgtitle('Переходные характеристики h(t)');
figure(2); sgtitle('Импульсные характеристики w(t)');

figure(3); clf; hold on; grid on;
T_values = T*[1/5 1/2 1 2 5];
fprintf('\nВлияние T на звено №4:\n');
for T_test = T_values
    W_test = K0/(T_test*s_tf + 1);
    plot(time, step(W_test, time), ...
        'DisplayName', ['T = ', num2str(T_test)]);
    t95 = -T_test*log(0.05);
    fprintf('T = %.2f: t95 = %.3f с\n', T_test, t95);
end
title('Влияние T на переходную характеристику звена №4');
xlabel('Время t, с'); ylabel('h(t)');
legend('show'); hold off;

figure(4); clf; hold on; grid on;
T1_values = T1*[1 1/2 1/5];
time_T1 = 0:0.01:40;
fprintf('\nВлияние T1 на звено №6:\n');
for T1_test = T1_values
    W_test = K0/(T2^2*s_tf^2 + T1_test*s_tf + 1);
    plot(time_T1, step(W_test, time_T1), ...
        'DisplayName', ['T1 = ', num2str(T1_test)]);
    info = stepinfo(W_test);
    tau = 2*T2^2/T1_test;
    fprintf('T1 = %.2f: перерегулирование = %.2f%%, tau = %.2f с\n', ...
        T1_test, info.Overshoot, tau);
end
title('Влияние T1 на колебательное звено №6 (T2 = const)');
xlabel('Время t, с'); ylabel('h(t)');
legend('show'); hold off;

fprintf('\nПроверка результатов при K = %.1f (максимальное расхождение):\n', K0);
for n = 1:11
    if ismember(n, [2 10])
        fprintf('Звено %2d: step/impulse не проверяются (идеальное дифференцирование)\n', n);
        continue;
    end

    [num, den] = numden(subs(W_all(n), K, K0));
    W_check = tf(sym2poly(num), sym2poly(den));

    y_step = step(W_check, time);
    h_num = double(subs(subs(h_all{n}, K, K0), t, time)) + zeros(size(time));
    err_step = max(abs(y_step(:) - h_num(:)));
    fprintf('Звено %2d: step = %.2e', n, err_step);

    if ismember(n, [1 4 5 6 7 8])
        y_impulse = impulse(W_check, time);
        w_num = double(subs(subs(w_all{n}, K, K0), t, time)) + zeros(size(time));
        err_impulse = max(abs(y_impulse(:) - w_num(:)));
        fprintf(', impulse = %.2e', err_impulse);
    else
        fprintf(', impulse пропущена (есть dirac-составляющая)');
    end
    fprintf('\n');
end
