clear;
clc;

K = 1;
T = 0.2;
T1 = 0.4;
T2 = 0.5;
T1_ap = 1.2;
K_values = K*[0.5 1 2];
time = (0:0.01:10).';
max_h_error = 0;
max_w_error = 0;

syms s t
Ts = sym(T, 'r');
T1s = sym(T1, 'r');
T2s = sym(T2, 'r');
T1aps = sym(T1_ap, 'r');

W_all = [
    1/s;
    s;
    1;
    1/(Ts*s + 1);
    1/(T2s^2*s^2 + T1aps*s + 1);
    1/(T2s^2*s^2 + T1s*s + 1);
    1/(Ts*s^2 + 1);
    1/(s*(Ts*s + 1));
    s/(Ts*s + 1);
    Ts*s + 1;
    (Ts*s + 1)/s
];

names = ["Идеальное интегрирующее", "Идеальное дифференцирующее", ...
    "Усилительное", "Апериодическое 1-го порядка", ...
    "Апериодическое 2-го порядка", "Колебательное", ...
    "Консервативное", "Инерционное интегрирующее", ...
    "Реальное дифференцирующее", "Форсирующее 1-го порядка", "Изодромное"];

U_step = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);

for j = 1:numel(W_all)
    W = W_all(j);
    Y_step = W*U_step;
    Y_impulse = W*U_impulse;

    D1 = limit(W/s, s, Inf);
    D0 = limit(W - D1*s, s, Inf);
    R = simplifyFraction(W - D1*s - D0);

    h = ilaplace(simplifyFraction(Y_step - D1), s, t);
    w = ilaplace(simplifyFraction(Y_impulse - D1*s - D0), s, t);
    h_values = valuesAtTime(h, t, time);
    w_values = valuesAtTime(w, t, time);

    [num, den] = numden(R);
    W_check = tf(K*sym2poly(num), sym2poly(den));
    h_check = step(W_check, time) + K*double(D0);
    w_check = impulse(W_check, time);
    max_h_error = max(max_h_error, max(abs(K*h_values - h_check(:))));
    max_w_error = max(max_w_error, max(abs(K*w_values - w_check(:))));

    figure('Name', names(j));
    subplot(1, 2, 1);
    hold on;
    grid on;
    for k = K_values
        plot(time, k*h_values, 'DisplayName', sprintf('K = %g', k));
    end
    title('Переходная h(t)');
    xlabel('Время t, с');
    legend('show', 'Location', 'best');
    hold off;

    subplot(1, 2, 2);
    hold on;
    grid on;
    for k = K_values
        plot(time, k*w_values, 'DisplayName', sprintf('K = %g', k));
    end
    title('Импульсная w(t)');
    xlabel('Время t, с');
    legend('show', 'Location', 'best');
    hold off;
    sgtitle(sprintf('%d. %s', j, names(j)));
end
fprintf('Максимальная ошибка проверки h(t): %.3e\n', max_h_error);
fprintf('Максимальная ошибка проверки w(t): %.3e\n', max_w_error);

time4 = 0:0.002:4;
figure('Name', 'Влияние T');
hold on;
grid on;
fprintf('\nЗвено 4. T, с | Время установления, с\n');
for Ti = T*[1/5 1/2 1 2 5]
    W_check = tf(K, [Ti 1]);
    y = step(W_check, time4);
    plot(time4, y, 'DisplayName', sprintf('T = %g', Ti));
    fprintf('%8.2f %15.3f\n', Ti, settlingTime(time4, y, K));
end
yline(0.95*K, ':', '95%', 'HandleVisibility', 'off');
title('Звено № 4. Влияние T');
xlabel('Время t, с');
legend('show', 'Location', 'best');

time6 = 0:0.01:50;
limits = [12 50];
T1_values = T1*[1 1/2 1/4 1/8 0];
figure('Name', 'Влияние T1');

y6 = zeros(numel(time6), numel(T1_values));
fprintf('\nЗвено 6. T1, с | Перерегулирование, %% | Время установления, с\n');
for j = 1:numel(T1_values)
    T1i = T1_values(j);
    W_check = tf(K, [T2^2 T1i 1]);
    y6(:,j) = step(W_check, time6);
    if T1i == 0
        fprintf('%8.2f %22s %22s\n', T1i, '—', 'Inf');
    else
        overshoot = 100*(max(y6(:,j)) - K)/K;
        t5 = settlingTime(time6, y6(:,j), K);
        fprintf('%8.2f %22.2f %22.3f\n', T1i, overshoot, t5);
    end
end

for panel = 1:2
    subplot(1, 2, panel);
    hold on;
    grid on;
    for j = 1:numel(T1_values)
        plot(time6, y6(:,j), 'DisplayName', sprintf('T1 = %g', T1_values(j)));
    end
    xlim([0 limits(panel)]);
    yline(0.95*K, ':', 'HandleVisibility', 'off');
    yline(1.05*K, ':', 'HandleVisibility', 'off');
    xlabel('Время t, с');
    legend('show', 'Location', 'best');
end
sgtitle('Звено № 6. Влияние T1 при фиксированном T2');

function t5 = settlingTime(time, y, finalValue)
    deviation = abs(y - finalValue)/abs(finalValue);
    last = find(deviation > 0.05, 1, 'last');
    if isempty(last)
        t5 = time(1);
    elseif last == numel(time)
        t5 = Inf;
    else
        t5 = interp1(deviation(last:last+1), time(last:last+1), 0.05);
    end
end

function y = valuesAtTime(formula, t, time)
    f = matlabFunction(formula, 'Vars', t);
    y = real(f(time)) + zeros(size(time));
    y(1) = double(limit(formula, t, 0, 'right'));
end
