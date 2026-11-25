clc;
clear;
close all;

s = tf('s');
K = 2.0; T = 0.5; T1 = 0.5; T2 = 0.4; T1a = 0.8;

names = {'1. Интегрирующее', '2. Дифференцирующее', '3. Усилительное', ...
    '4. Апериодическое 1-го порядка', '5. Апериодическое 2-го порядка', ...
    '6. Колебательное', '7. Консервативное', '8. Инерционное интегрирующее', ...
    '9. Реальное дифференцирующее', '10. Форсирующее', '11. Изодромное'};

W = {K/s, K*s, tf(K), K/(T*s + 1), K/(T2^2*s^2 + T1a*s + 1), ...
    K/(T2^2*s^2 + T1*s + 1), K/(T*s^2 + 1), K/(s*(T*s + 1)), ...
    K*s/(T*s + 1), K*(T*s + 1), K*(T*s + 1)/s};

w = logspace(-2, 2, 2000);
n = numel(names);

for i = 1:n
    figure;
    nyquist(W{i});
    grid on;
    title(names{i});

    figure;
    bode(W{i});
    grid on;
    title(names{i});
end

K_values = [0.5 1 2 4];
for i = 1:n
    figure;
    for Kv = K_values
        sys = Kv*(W{i}/K);
        [mag, ph] = bode(sys, w);
        subplot(2, 1, 1);
        semilogx(w, 20*log10(squeeze(mag)));
        hold on; grid on;
        subplot(2, 1, 2);
        semilogx(w, squeeze(ph));
        hold on; grid on;
    end
end

T_values = [0.1 0.25 0.5 1 2.5];
T2_values = [0.3 0.4 0.6 0.8];
for i = 1:n
    figure;
    switch i
        case 1
            models = {K/s};
        case 2
            models = {K*s};
        case 3
            models = {tf(K)};
        case 4
            models = arrayfun(@(x) K/(x*s + 1), T_values, 'UniformOutput', false);
        case 5
            models = arrayfun(@(x) K/(x^2*s^2 + T1a*s + 1), T2_values, 'UniformOutput', false);
        case 6
            models = arrayfun(@(x) K/(x^2*s^2 + T1*s + 1), T2_values, 'UniformOutput', false);
        case 7
            models = arrayfun(@(x) K/(x*s^2 + 1), T_values, 'UniformOutput', false);
        case 8
            models = arrayfun(@(x) K/(s*(x*s + 1)), T_values, 'UniformOutput', false);
        case 9
            models = arrayfun(@(x) K*s/(x*s + 1), T_values, 'UniformOutput', false);
        case 10
            models = arrayfun(@(x) K*(x*s + 1), T_values, 'UniformOutput', false);
        case 11
            models = arrayfun(@(x) K*(x*s + 1)/s, T_values, 'UniformOutput', false);
    end
    for j = 1:numel(models)
        [mag, ph] = bode(models{j}, w);
        subplot(2, 1, 1);
        semilogx(w, 20*log10(squeeze(mag)));
        hold on; grid on;
        subplot(2, 1, 2);
        semilogx(w, squeeze(ph));
        hold on; grid on;
    end
end
