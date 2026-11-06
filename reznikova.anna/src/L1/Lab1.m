syms s t K;
T = 1.5; T1 = 1.6; T2 = 2.0;
T1a = 5;

names = {'интегрирующего звена', 'дифференцирующего звена', 'усилительного звена', ...
    'апериодического звена 1-го порядка', 'апериодического звена 2-го порядка', ...
    'колебательного звена', 'консервативного звена', 'инерционного интегрирующего звена', ...
    'реального дифференцирующего звена', 'форсирующего звена', 'изодромного звена'};

n = numel(names);

W = {K/s, K*s, K, K/(T*s+1), K/(T2^2*s^2+T1a*s+1), K/(T2^2*s^2+T1*s+1), ...
    K/(T*s^2+1), K/(s*(T*s+1)), K*s/(T*s+1), K*(T*s+1), K*(T*s+1)/s};

U_step = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);

h = cell(1, n);
w = cell(1, n);

for i = 1:n
    Y_step = W{i}*U_step;
    h{i} = ilaplace(Y_step, s, t);
    Y_impulse = W{i}*U_impulse;
    w{i} = ilaplace(Y_impulse, s, t);

    fprintf('%d. %s\n', i, names{i});
    disp(h{i});
    disp(w{i});
end

K_values = 4.0;
time = 0:0.01:20;

s_tf = tf('s');
W_tf = {@(k) k/s_tf, [], @(k) tf(k, 1), @(k) k/(T*s_tf+1), ...
    @(k) k/(T2^2*s_tf^2+T1a*s_tf+1), @(k) k/(T2^2*s_tf^2+T1*s_tf+1), ...
    @(k) k/(T*s_tf^2+1), @(k) k/(s_tf*(T*s_tf+1)), @(k) k*s_tf/(T*s_tf+1), ...
    [], @(k) k*(T*s_tf+1)/s_tf};

for i = 1:n
    figure('Position', [100 100 1100 700]);

    subplot(2, 2, 1);
    hold on;
    grid on;
    for k = K_values
        h_k = subs(h{i}, K, k);
        y = double(subs(h_k, t, time)) + zeros(size(time));
        plot(time, y, 'DisplayName', ['K = ', num2str(k)]);
    end
    title(['Переходные характеристики ', names{i}]);
    xlabel('Время t, с');
    ylabel('h(t)');
    legend('show', 'Location', 'northwest');
    hold off;

    subplot(2, 2, 2);
    hold on;
    grid on;
    for k = K_values
        w_k = subs(w{i}, K, k);
        y = double(subs(w_k, t, time)) + zeros(size(time));
        plot(time, y, 'DisplayName', ['K = ', num2str(k)]);
    end
    title(['Импульсные характеристики ', names{i}]);
    xlabel('Время t, с');
    ylabel('w(t)');
    legend('show', 'Location', 'northwest');
    hold off;

    if isempty(W_tf{i})
        continue
    end

    subplot(2, 2, 3);
    hold on;
    grid on;
    for k = K_values
        W_check = W_tf{i}(k);
        y = step(W_check, time);
        plot(time, y, 'DisplayName', ['K = ', num2str(k)]);
    end
    title('Проверка переходной характеристики');
    xlabel('Время t, с');
    ylabel('Выходной сигнал');
    legend('show', 'Location', 'northwest');
    hold off;

    subplot(2, 2, 4);
    hold on;
    grid on;
    for k = K_values
        W_check = W_tf{i}(k);
        y = impulse(W_check, time);
        plot(time, y, 'DisplayName', ['K = ', num2str(k)]);
    end
    title('Проверка импульсной характеристики');
    xlabel('Время t, с');
    ylabel('Выходной сигнал');
    legend('show', 'Location', 'northwest');
    hold off;
end

K0 = 4.0;
time_star = 0:0.01:40;

figure;
hold on;
grid on;
for Ti = T*[0.2 0.5 1 2 5]
    W_star = tf(K0, [Ti 1]);
    y = step(W_star, time_star);
    plot(time_star, y, 'DisplayName', ['T = ', num2str(Ti)]);
    info = stepinfo(W_star);
    fprintf('T = %.2f, время переходного процесса = %.2f с\n', Ti, info.SettlingTime);
end
title('Влияние T на переходную характеристику апериодического звена');
xlabel('Время t, с');
ylabel('h(t)');
legend('show', 'Location', 'southeast');
hold off;

time_star = 0:0.01:60;

figure;
hold on;
grid on;
for T1i = [1.6 0.8 0.4 0.2 0.1]
    W_star = tf(K0, [T2^2 T1i 1]);
    y = step(W_star, time_star);
    plot(time_star, y, 'DisplayName', ['T1 = ', num2str(T1i)]);
    info = stepinfo(W_star);
    fprintf('T1 = %.2f, перерегулирование = %.1f %%, время затухания = %.1f с\n', ...
        T1i, info.Overshoot, info.SettlingTime);
end
title('Влияние T1 на переходную характеристику колебательного звена');
xlabel('Время t, с');
ylabel('h(t)');
legend('show', 'Location', 'northeast');
hold off;
