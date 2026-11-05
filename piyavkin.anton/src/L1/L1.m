clear; clc; close all;

K   = 2.8;
T   = 0.6;
T1  = 0.6;
T2  = 0.5;
T1a = 1.5;

t_end = 20;
dt    = 0.005;
t     = 0:dt:t_end;

link_names = {
    '1. Интегрирующее'
    '2. Дифференцирующее'
    '3. Усилительное'
    '4. Апериодическое 1-го порядка'
    '5. Апериодическое 2-го порядка'
    '6. Колебательное'
    '7. Консервативное'
    '8. Инерционное интегрирующее'
    '9. Реальное дифференцирующее'
    '10. Форсирующее 1-го порядка'
    '11. Изодромное'
};

syms ss tau k real

W_list = {
    k/ss
    k*ss
    k
    k/(T*ss + 1)
    k/(T2^2*ss^2 + T1a*ss + 1)
    k/(T2^2*ss^2 + T1*ss + 1)
    k/(T*ss^2 + 1)
    k/(ss*(T*ss + 1))
    k*ss/(T*ss + 1)
    k*(T*ss + 1)
    k*(T*ss + 1)/ss
};

Ustep = laplace(heaviside(tau), tau, ss);
Uimp = laplace(dirac(tau), tau, ss);

fprintf('Аналитические выражения временных характеристик\n');
fprintf('K = %.1f, T = %.2f, T1 = %.2f, T2 = %.2f, T1a = %.2f\n\n', K, T, T1, T2, T1a);

for n = 1:numel(W_list)
    Ws = subs(W_list{n}, k, K);
    Hs = simplify(ilaplace(Ws * Ustep, ss, tau));
    Ws_t = simplify(ilaplace(Ws * Uimp,  ss, tau));
    fprintf('%s\n', link_names{n});
    fprintf('  h(t) = %s\n', char(Hs));
    fprintf('  w(t) = %s\n\n', char(Ws_t));
end

s = tf('s');
sys = {
    K/s
    K*s
    tf(K)
    K/(T*s + 1)
    K/(T2^2*s^2 + T1a*s + 1)
    K/(T2^2*s^2 + T1*s + 1)
    K/(T*s^2 + 1)
    K/(s*(T*s + 1))
    K*s/(T*s + 1)
    K*(T*s + 1)
    K*(T*s + 1)/s
};

figure('Name','h(t)','Color','w','Position',[80 80 1100 720]);
for n = 1:numel(sys)
    subplot(3,4,n); grid on; hold on; box on;
    switch n
        case 2
            plot(t, zeros(size(t)), 'b', 'LineWidth', 1.3);
            stem(0, K, 'b', 'filled');
        case 10
            plot(t, K*ones(size(t)), 'b', 'LineWidth', 1.3);
            stem(0, K*T, 'b', 'filled');
        otherwise
            [y, ty] = step(sys{n}, t);
            plot(ty, y, 'b', 'LineWidth', 1.3);
    end
    title(link_names{n}, 'FontSize', 8);
    xlabel('t, с'); ylabel('h(t)');
end
sgtitle(sprintf('Переходные функции, K = %.1f', K));

col_w = [0.75 0.1 0.1];
figure('Name','w(t)','Color','w','Position',[100 60 1100 720]);
for n = 1:numel(sys)
    subplot(3,4,n); grid on; hold on; box on;
    switch n
        case 2
            plot(t, zeros(size(t)), 'Color',col_w, 'LineWidth', 1.3);
            stem([-dt dt], [K -K], 'Color',col_w, 'filled');
        case 3
            plot(t, zeros(size(t)), 'Color',col_w, 'LineWidth', 1.3);
            stem(0, K, 'Color',col_w, 'filled');
        case 9
            yr = -(K/T^2)*exp(-t/T);
            plot(t, yr, 'Color',col_w, 'LineWidth', 1.3);
            stem(0, K/T, 'Color',col_w, 'filled');
        case 10
            plot(t, zeros(size(t)), 'Color',col_w, 'LineWidth', 1.3);
            stem([-dt 0 dt], [K*T  K  -K*T], 'Color',col_w, 'filled');
        case 11
            plot(t, K*ones(size(t)), 'Color',col_w, 'LineWidth', 1.3);
            stem(0, K*T, 'Color',col_w, 'filled');
        otherwise
            [y, ty] = impulse(sys{n}, t);
            plot(ty, y, 'Color',col_w, 'LineWidth', 1.3);
    end
    title(link_names{n}, 'FontSize', 8);
    xlabel('t, с'); ylabel('w(t)');
end
sgtitle(sprintf('Функции веса, K = %.1f', K));

figure('Name','Влияние K','Color','w'); hold on; grid on; box on;
Kset = [1.4 2.8 5.6];
cols  = lines(numel(Kset));
for i = 1:numel(Kset)
    G = Kset(i)/(T2^2*s^2 + T1*s + 1);
    [y, ty] = step(G, t);
    plot(ty, y, 'Color', cols(i,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('K = %.1f', Kset(i)));
end
title('Колебательное звено: изменение K');
xlabel('t, с'); ylabel('h(t)');
legend('Location','best'); legend boxoff;

figure('Name','Влияние T на h(t)','Color','w'); hold on; grid on; box on;
Tset = T*[0.2 0.5 1 2 5];
t2   = 0:0.01:max(8*max(Tset), 15);
for i = 1:numel(Tset)
    G = K/(Tset(i)*s + 1);
    [y, ty] = step(G, t2);
    plot(ty, y, 'LineWidth', 1.5, 'DisplayName', sprintf('T = %.2g с', Tset(i)));
end
title('Апериодическое звено №4: влияние T');
xlabel('t, с'); ylabel('h(t)');
legend('Location','best'); legend boxoff;

figure('Name','Влияние T на w(t)','Color','w'); hold on; grid on; box on;
for i = 1:numel(Tset)
    G = K/(Tset(i)*s + 1);
    [y, ty] = impulse(G, t2);
    plot(ty, y, 'LineWidth', 1.5, 'DisplayName', sprintf('T = %.2g с', Tset(i)));
end
title('Апериодическое звено №4: импульсные характеристики');
xlabel('t, с'); ylabel('w(t)');
legend('Location','best'); legend boxoff;

figure('Name','Влияние T1','Color','w'); hold on; grid on; box on;
T1set = [0.60 0.40 0.20 0.10 0.04];
t3    = 0:0.01:30;
for i = 1:numel(T1set)
    G = K/(T2^2*s^2 + T1set(i)*s + 1);
    [y, ty] = step(G, t3);
    plot(ty, y, 'LineWidth', 1.5, 'DisplayName', sprintf('T_1 = %.2g с', T1set(i)));
end
title('Колебательное звено: уменьшение T_1 при T_2 = 0.5 с');
xlabel('t, с'); ylabel('h(t)');
legend('Location','best'); legend boxoff;
xlim([0 20]);

fprintf('Расчёт завершён. Открыты окна с графиками.\n');
