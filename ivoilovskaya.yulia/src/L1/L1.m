clear; clc; close all;
syms s t

K = 2.2;
T = 0.2;
T1 = 0.4;
T2 = 0.3;
K_values = [K/2 K 2*K];
time = 0:0.01:5;
colors = [0 0.35 0.70; 0.85 0.33 0.10; 0.15 0.55 0.25];

names = ["Интегрирующее", "Идеальное дифференцирующее", ...
    "Усилительное", "Апериодическое 1-го порядка", ...
    "Звено 2-го порядка, строка 5", "Колебательное", ...
    "Консервативное", "Инерционное интегрирующее", ...
    "Реальное дифференцирующее", "Форсирующее", "Изодромное"];

W(1) = K/s;
W(2) = K*s;
W(3) = K;
W(4) = K/(T*s + 1);
W(5) = K/(T2^2*s^2 + T1*s + 1);
W(6) = K/(T2^2*s^2 + T1*s + 1);
W(7) = K/(T*s^2 + 1);
W(8) = K/(s*(T*s + 1));
W(9) = K*s/(T*s + 1);
W(10) = K*(T*s + 1);
W(11) = K*(T*s + 1)/s;

p = tf('s');
models{1} = K/p;
models{2} = K*p;
models{3} = tf(K);
models{4} = K/(T*p + 1);
models{5} = K/(T2^2*p^2 + T1*p + 1);
models{6} = K/(T2^2*p^2 + T1*p + 1);
models{7} = K/(T*p^2 + 1);
models{8} = K/(p*(T*p + 1));
models{9} = K*p/(T*p + 1);
models{10} = K*(T*p + 1);
models{11} = K*(T*p + 1)/p;

U_step = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);
h_note = strings(1,11);
w_note = strings(1,11);
h_note(2) = "Обычная часть равна нулю. Полный отклик: K δ(t).";
h_note(10) = "К обычной части добавить K T δ(t).";
w_note(2) = "Обычная часть равна нулю. Полный отклик: K δ'(t).";
w_note(3) = "Обычная часть равна нулю. Полный отклик: K δ(t).";
w_note(9) = "К обычной части добавить (K/T) δ(t).";
w_note(10) = "Обычная часть равна нулю. Отклик: K δ(t) + K T δ'(t).";
w_note(11) = "К обычной части добавить K T δ(t).";
disp('При T1 = 0.4 < 2*T2 = 0.6 звенья 5 и 6 совпадают.');

for n = 1:11
    h = ilaplace(W(n)*U_step, s, t);
    w = ilaplace(W(n)*U_impulse, s, t);
    disp(names(n));
    disp('h(t) при K = 2.2:'); disp(h);
    disp('w(t) при K = 2.2:'); disp(w);

    h_plot = subs(h, dirac(t), 0);
    w_plot = subs(w, dirac(1,t), 0);
    w_plot = subs(w_plot, dirac(t), 0);
    h_values = double(subs(h_plot, t, time));
    w_values = double(subs(w_plot, t, time));
    if isscalar(h_values), h_values = h_values*ones(size(time)); end
    if isscalar(w_values), w_values = w_values*ones(size(time)); end

    if n ~= 2 && n ~= 10
        h_check = step(models{n}, time);
        w_check = impulse(models{n}, time);
        fprintf('Проверка при K = 2.2: ошибка h = %.2e, w = %.2e\n', ...
            max(abs(h_check(:)-h_values(:))), max(abs(w_check(:)-w_values(:))));
    end

    figure('Name',names(n),'Color','w','DefaultTextInterpreter','none', ...
        'DefaultLegendInterpreter','none','DefaultAxesTickLabelInterpreter','none');
    theme(gcf, 'light');
    subplot(2,1,1); hold on; grid on;
    for j = 1:3
        plot(time, h_values*K_values(j)/K, 'Color', colors(j,:), ...
            'LineWidth', 1.8, 'DisplayName', sprintf('K = %.1f', K_values(j)));
    end
    title('Переходная характеристика (обычная часть)');
    xlabel('t, с'); ylabel('h(t)'); xlim([0 5]);
    legend('show', 'Location', 'eastoutside');
    if strlength(h_note(n)) > 0
        subtitle(h_note(n),'Interpreter','none');
    end

    subplot(2,1,2); hold on; grid on;
    for j = 1:3
        plot(time, w_values*K_values(j)/K, 'Color', colors(j,:), ...
            'LineWidth', 1.8, 'DisplayName', sprintf('K = %.1f', K_values(j)));
    end
    title('Импульсная характеристика (обычная часть)');
    xlabel('t, с'); ylabel('w(t)'); xlim([0 5]);
    legend('show', 'Location', 'eastoutside');
    if strlength(w_note(n)) > 0
        subtitle(w_note(n),'Interpreter','none');
    end
    sgtitle(sprintf('%d. %s', n, names(n)));
end

T_values = [T/5 T/2 T 2*T 5*T];
T_colors = [colors; 0.49 0.18 0.56; 0.2 0.2 0.2];
figure('Name','Влияние T','Color','w','DefaultTextInterpreter','none', ...
    'DefaultLegendInterpreter','none','DefaultAxesTickLabelInterpreter','none');
theme(gcf,'light');
hold on; grid on;
for j = 1:5
    currentT = T_values(j);
    model = tf(K, [currentT 1]);
    plot(time, step(model,time), 'Color',T_colors(j,:),'LineWidth',1.8, ...
        'DisplayName', sprintf('T = %.2f с',currentT));
    info = stepinfo(model);
    fprintf('T = %.2f: время установления = %.3f с\n',currentT,info.SettlingTime);
end
title('Апериодическое звено: влияние T при K = 2.2');
xlabel('t, с'); ylabel('h(t)'); xlim([0 5]);
legend('show','Location','southeast');

T1_values = [T1 T1/2 T1/5];
time_long = 0:0.01:15;
figure('Name','Влияние T1','Color','w','DefaultTextInterpreter','none', ...
    'DefaultLegendInterpreter','none','DefaultAxesTickLabelInterpreter','none');
theme(gcf,'light');
subplot(2,1,1); hold on; grid on;
for j = 1:3
    currentT1 = T1_values(j);
    model = tf(K, [T2^2 currentT1 1]);
    plot(time_long, step(model,time_long), 'Color',colors(j,:), ...
        'LineWidth',1.8,'DisplayName',sprintf('T1 = %.2f с',currentT1));
    info = stepinfo(model);
    fprintf('T1 = %.2f: перерегулирование = %.2f%%, время установления = %.3f с\n', ...
        currentT1,info.Overshoot,info.SettlingTime);
end
title('Уменьшение T1: затухающие колебания');
xlabel('t, с'); ylabel('h(t)'); xlim([0 15]);
legend('show','Location','northeast');

subplot(2,1,2);
time_zero = 0:0.01:6;
model = tf(K, [T2^2 0 1]);
plot(time_zero,step(model,time_zero),'Color',colors(1,:),'LineWidth',1.8);
grid on;
title('T1 = 0: незатухающие колебания');
xlabel('t, с'); ylabel('h(t)'); xlim([0 6]);
sgtitle('Колебательное звено: T2 = 0.3 с, K = 2.2');
disp('При T1 = 0 колебания не затухают, время установления бесконечно.');
