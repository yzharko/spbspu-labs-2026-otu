syms s t K;
K0 = 1.2;
T = 0.3;
T1 = 0.5;
T2 = 0.6;
K_values = K0*[1 2 3 4];
time = 0:0.01:10;
s_tf = tf('s');

names = {'Интегрирующее', 'Идеальное дифференцирующее', ...
    'Усилительное', 'Апериодическое 1-го порядка', ...
    'Звено 5: условие T1 >= 2*T2', 'Колебательное', ...
    'Консервативное', 'Инерционное интегрирующее', ...
    'Реальное дифференцирующее', 'Форсирующее', 'Изодромное'};

W_all = [K/s; K*s; K; K/(T*s+1); ...
    K/(T2^2*s^2+T1*s+1); K/(T2^2*s^2+T1*s+1); ...
    K/(T*s^2+1); K/(s*(T*s+1)); K*s/(T*s+1); ...
    K*(T*s+1); K*(T*s+1)/s];


window_names = {'Переходные характеристики: регулярные части', ...
    'Импульсные характеристики: регулярные части'};
windows = cell(1, 2);
labels = {'h(t)', 'w(t)'};
colors = lines(numel(K_values));
for p = 1:2
    windows{p} = figure('Name', window_names{p}, ...
        'NumberTitle', 'off', 'Position', [50, 50, 1400, 850]);
end


U_step = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);

for n = 1:11
    W = W_all(n);
    h = ilaplace(W*U_step, s, t);
    w = ilaplace(W*U_impulse, s, t);

    h_plot = subs(h, dirac(t), 0);
    w_plot = subs(w, dirac(1, t), 0);
    w_plot = subs(w_plot, dirac(t), 0);

    ax = cell(1, 2);
    for p = 1:2
        figure(windows{p});
        ax{p} = subplot(3, 4, n);
        hold(ax{p}, 'on'); 
        grid(ax{p}, 'on');
        
        xlabel(ax{p}, 't, с'); 
        ylabel(ax{p}, labels{p});
    end

    for j = 1:numel(K_values)
        k = K_values(j);
        h_k = subs(h_plot, K, k);
        w_k = subs(w_plot, K, k);
        y_h = double(subs(h_k, t, time)) + zeros(size(time));
        y_w = double(subs(w_k, t, time)) + zeros(size(time));

        plot(ax{1}, time, y_h, 'Color', colors(j,:));
        plot(ax{2}, time, y_w, 'Color', colors(j,:));

    end

    title(ax{1}, {sprintf('%d. %s', n, names{n}), ...
        ['Добавка: ' char(simplify(h-h_plot))]}, ...
        'Interpreter', 'none', 'FontSize', 8);
    title(ax{2}, {sprintf('%d. %s', n, names{n}), ...
        ['Добавка: ' char(simplify(w-w_plot))]}, ...
        'Interpreter', 'none', 'FontSize', 8);
end

for p = 1:2
    figure(windows{p});
    subplot(3, 4, 12); hold on;
    for j = 1:numel(K_values)
        plot(NaN, NaN, 'Color', colors(j,:), ...
            'DisplayName', sprintf('K = %g', K_values(j)));
    end
    axis off;
    title('Значения K');
    legend('show', 'Location', 'northwest', 'FontSize', 8);
end

figure('Name', 'Влияние T: апериодическое звено 1-го порядка');
hold on;
grid on;
T_values = T*[1/5, 1/2, 1, 2, 5];
for T_test = T_values
    W_test = K0/(T_test*s_tf+1);
    plot(time, step(W_test, time), ...
        'DisplayName', sprintf('T = %g', T_test));
end

title('Влияние T на h(t)');
xlabel('t, с'); 
ylabel('h(t)'); 
legend('show');
hold off;

figure('Name', 'Влияние T1: колебательное звено');
hold on;
grid on;
T1_values = T1*[1, 1/2, 1/5];
time_T1 = 0:0.01:40; 
for T1_test = T1_values
    W_test = K0/(T2^2*s_tf^2+T1_test*s_tf+1);
    plot(time_T1, step(W_test, time_T1), ...
        'DisplayName', sprintf('T1 = %g', T1_test));
end

title('Влияние T1 на h(t), T2=const');
xlabel('t, с'); 
ylabel('h(t)'); 
legend('show');
hold off;
