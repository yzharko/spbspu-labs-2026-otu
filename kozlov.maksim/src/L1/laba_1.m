clear; clc; close all;

K = 4.0; 
T = 1.5; 
T1 = 1.6; 
T2 = 2.0;
T1a = 4.5;
dt = 0.01;

time = 0:dt:30; 
s = tf('s');

W = cell(11,1);
W{1} = K/s;
W{2} = K*s;
W{3} = K;
W{4} = K/(T*s + 1);
W{5} = K/(T2^2*s^2 + T1a*s + 1);
W{6} = K/(T2^2*s^2 + T1*s + 1);
W{7} = K/(T*s^2 + 1);
W{8} = K/(s*(T*s + 1));
W{9} = K*s/(T*s + 1);
W{10} = K*(T*s + 1);
W{11} = K*(T*s + 1)/s;

names = {'1. Интегрирующее', '2. Дифференцирующее', '3. Усилительное', ...
         '4. Апериодическое 1-го', '5. Апериодическое 2-го', '6. Колебательное', ...
         '7. Консервативное', '8. Инерционное интегр.', '9. Реальное дифф.', ...
         '10. Форсирующее', '11. Изодромное'};

s_sym = sym('s');
t_sym = sym('t');

W_sym = cell(11,1);
W_sym{1} = K/s_sym;
W_sym{2} = K*s_sym;
W_sym{3} = K;
W_sym{4} = K/(T*s_sym + 1);
W_sym{5} = K/(T2^2*s_sym^2 + T1a*s_sym + 1);
W_sym{6} = K/(T2^2*s_sym^2 + T1*s_sym + 1);
W_sym{7} = K/(T*s_sym^2 + 1);
W_sym{8} = K/(s_sym*(T*s_sym + 1));
W_sym{9} = K*s_sym/(T*s_sym + 1);
W_sym{10} = K*(T*s_sym + 1);
W_sym{11} = K*(T*s_sym + 1)/s_sym;

fig1 = figure('Name', 'h(t)', 'Color', 'w', 'Position', [50 50 1400 900]);
for i = 1:numel(names)
    subplot(3,4,i); hold on; grid on;
    
    if i == 2
        y = zeros(size(time)); 
        y(1) = K / dt; 
        plot(time, y, 'b');
    elseif i == 3
        y = K * ones(size(time));
        plot(time, y, 'b');
    elseif i == 10
        y = K * ones(size(time));
        y(1) = K + (K * T) / dt;
        plot(time, y, 'b');
    else
        [y, t_out] = step(W{i}, time);
        plot(t_out, y, 'b');
    end
    
    Y_s = W_sym{i} * (1/s_sym);
    y_t = ilaplace(Y_s, s_sym, t_sym);
    y_analyt = double(subs(y_t, t_sym, time));
    plot(time, y_analyt, 'r--', 'LineWidth', 1.5);
    
    xlim([0 30]);
    title(names{i}, 'FontSize', 11, 'FontWeight', 'bold');
    xlabel('t, c', 'FontSize', 11);
    ylabel('h(t)', 'FontSize', 11);
end
exportgraphics(fig1, '01_Переходные_h.png');

fig2 = figure('Name', 'w(t)', 'Color', 'w', 'Position', [100 100 1400 900]);
for i = 1:numel(names)
    subplot(3,4,i); hold on; grid on;
    
    if i == 2
        y = zeros(size(time));
        y(1) = K / (dt^2); 
        y(2) = -K / (dt^2);
        plot(time, y, 'r');
    elseif i == 3
        y = zeros(size(time)); 
        y(1) = K / dt;
        plot(time, y, 'r');
    elseif i == 10
        y = zeros(size(time)); 
        y(1) = (K * T) / (dt^2); 
        y(2) = -(K * T) / (dt^2);
        plot(time, y, 'r');
    elseif i == 11
        y = K * ones(size(time)); 
        y(1) = K + (K * T) / dt;
        plot(time, y, 'r');
        ylim([0 10]);
    else
        [y, t_out] = impulse(W{i}, time);
        plot(t_out, y, 'r');
    end
    
    Y_s = W_sym{i};
    y_t = ilaplace(Y_s, s_sym, t_sym);
    y_analyt = double(subs(y_t, t_sym, time));
    plot(time, y_analyt, 'm--', 'LineWidth', 1.5);
    
    xlim([0 30]); 
    title(names{i}, 'FontSize', 11, 'FontWeight', 'bold');
    xlabel('t, c', 'FontSize', 11);
    ylabel('w(t)', 'FontSize', 11);
end
exportgraphics(fig2, '02_Импульсные_w.png');

fig3 = figure('Name', 'Влияние K', 'Color', 'w', 'Position', [150 150 900 600]);
hold on; grid on;
K_vals = [1.0, 2.0, 4.0, 6.0, 8.0];
for k = K_vals
    [y, t_out] = step(k/(T2^2*s^2 + T1*s + 1), time);
    plot(t_out, y); 
end
title('Колебательное звено: варьирование K', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('t, c', 'FontSize', 12); ylabel('h(t)', 'FontSize', 12);
legend('K=1','K=2','K=4','K=6','K=8', 'Location', 'northeast');
xlim([0 30]);
exportgraphics(fig3, '03_Влияние_K.png');

fig4 = figure('Name', 'Влияние T', 'Color', 'w', 'Position', [200 200 900 600]);
hold on; grid on;
T_factors = [0.2, 0.5, 1.0, 2.0, 5.0];
for f = T_factors
    [y, t_out] = step(K/((T*f)*s + 1), time);
    plot(t_out, y); 
end
title('Апериодическое звено: варьирование T', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('t, c', 'FontSize', 12); ylabel('h(t)', 'FontSize', 12);
legend('T/5','T/2','T','2T','5T', 'Location', 'southeast');
xlim([0 30]);
exportgraphics(fig4, '04_Влияние_T.png');

fig5 = figure('Name', 'Влияние T1', 'Color', 'w', 'Position', [250 250 900 600]);
hold on; grid on;
T1_vals = [1.6, 1.2, 0.8, 0.4, 0.1];
for t1 = T1_vals
    [y, t_out] = step(K/(T2^2*s^2 + t1*s + 1), time);
    plot(t_out, y); 
end
title('Колебательное звено: уменьшение T1', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('t, c', 'FontSize', 12); ylabel('h(t)', 'FontSize', 12);
legend('T1=1.6','T1=1.2','T1=0.8','T1=0.4','T1=0.1', 'Location', 'northeast');
xlim([0 30]);
exportgraphics(fig5, '05_Влияние_T1.png');

