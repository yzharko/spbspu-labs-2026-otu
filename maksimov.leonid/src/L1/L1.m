clear; clc; close all;

K0 = 2.4;
T0 = 0.4;
T1 = 0.8;
T2 = 0.9;

K_values = K0 * (1:4);
time = 0:0.01:10;
tt = time;  tt(1) = 1e-9;

syms s K T positive;
syms t real;

T1r = sym(T1, 'r');  T2r = sym(T2, 'r');

names = {'Интегрирующее', 'Дифференцирующее', 'Усилительное', ...
         'Апериодическое 1-го порядка', 'Апериодическое 2-го порядка', ...
         'Колебательное', 'Консервативное', 'Инерционное интегрирующее', ...
         'Реальное дифференцирующее', 'Форсирующее 1-го порядка', 'Изодромное'};

Wsym = { K/s, ...
         K*s, ...
         K, ...
         K/(T*s + 1), ...
         K/(T2r^2*s^2 + T1r*s + 1), ...
         K/(T2r^2*s^2 + T1r*s + 1), ...
         K/(T*s^2 + 1), ...
         K/(s*(T*s + 1)), ...
         K*s/(T*s + 1), ...
         K*(T*s + 1), ...
         K*(T*s + 1)/s };

U_step    = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);

for n = 1:numel(names)
    W = Wsym{n};

    Y_step = simplify(W*U_step);
    h = ilaplace(Y_step, s, t);

    Y_impulse = simplify(W*U_impulse);
    w = ilaplace(Y_impulse, s, t);

    fprintf('\n%d. %s\n', n, names{n});
    fprintf('W(s) = %s\n', char(W));
    fprintf('h(t) = %s\n', char(h));
    fprintf('w(t) = %s\n', char(w));

    hf = matlabFunction(h, 'Vars', [K T t]);
    wf = matlabFunction(w, 'Vars', [K T t]);

    figure('Name', names{n}, 'Position', [100 100 1100 450]);

    subplot(1, 2, 1); hold on; grid on;
    for k = K_values
        h_k = hf(k, T0, tt);
        y = real(h_k + zeros(size(tt)));
        plot(time, y, 'DisplayName', ['K = ', num2str(k)]);
    end
    title(['Переходные характеристики: ', names{n}]);
    xlabel('Время t, с');
    ylabel('h(t)');
    legend('show', 'Location', 'best');
    hold off;

    subplot(1, 2, 2); hold on; grid on;
    for k = K_values
        w_k = wf(k, T0, tt);
        y = real(w_k + zeros(size(tt)));
        plot(time, y, 'DisplayName', ['K = ', num2str(k)]);
    end
    title(['Импульсные характеристики: ', names{n}]);
    xlabel('Время t, с');
    ylabel('w(t)');
    legend('show', 'Location', 'best');
    hold off;
end

T_list = T0 * [1/5, 1/2, 1, 2, 5];
timeL = 0:0.01:20;

figure; hold on; grid on;
fprintf('\n--- Апериодическое звено (K = %g) ---\n', K0);
for Tc = T_list
    y = step(tf(K0, [Tc 1]), timeL);
    ts = settle_time(timeL, y, K0, 0.05);
    plot(timeL, y, 'DisplayName', ['T = ', num2str(Tc)]);
    fprintf('T = %-5g  время регулирования (5%%) = %.3f с\n', Tc, ts);
end
title('Апериодическое звено: влияние постоянной времени T');
xlabel('Время t, с');
ylabel('h(t)');
legend('show', 'Location', 'southeast');
hold off;

T1_list = T1 * [1, 1/2, 1/4, 1/8];
timeO = 0:0.01:100;

figure; hold on; grid on;
fprintf('\n--- Колебательное звено (K = %g, T2 = %g) ---\n', K0, T2);
for Tc1 = T1_list
    y = step(tf(K0, [T2^2 Tc1 1]), timeO);
    sigma = (max(y) - K0) / K0 * 100;
    ts = settle_time(timeO, y, K0, 0.05);
    plot(timeO, y, 'DisplayName', ['T1 = ', num2str(Tc1)]);
    fprintf('T1 = %-5g  перерегулирование = %.2f %%,  время затухания (5%%) = %.2f с\n', Tc1, sigma, ts);
end
title('Колебательное звено: влияние T1 при T2 = const');
xlabel('Время t, с');
ylabel('h(t)');
xlim([0 60]);
legend('show', 'Location', 'northeast');
hold off;

function ts = settle_time(t, y, yss, tol)
    idx = find(abs(y - yss) > tol*abs(yss), 1, 'last');
    if isempty(idx)
        ts = 0;
    elseif idx == numel(t)
        ts = NaN;
    else
        ts = t(idx + 1);
    end
end
