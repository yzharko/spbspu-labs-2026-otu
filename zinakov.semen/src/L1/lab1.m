clear; clc; close all;

K = 2.2; 
T = 0.2; 
T1 = 0.4; 
T2 = 0.3;

% для апериодического звена 2-го порядка, тк там надо T1 >= 2*T2 = 0.6 с
T1a = 0.7; 

dt = 0.001; 
time = 0:dt:5; 
s = tf('s');

names = {'1. Интегрирующее', '2. Дифференцирующее', '3. Усилительное', ...
         '4. Апериодическое 1-го', '5. Апериодическое 2-го', '6. Колебательное', ...
         '7. Консервативное', '8. Инерционное интегр.', '9. Реальное дифф.', ...
         '10. Форсирующее', '11. Изодромное'};

W = cell(11,1);
W{1}  = K/s;
W{2}  = K*s;
W{3}  = K;
W{4}  = K/(T*s + 1);
W{5}  = K/(T2^2*s^2 + T1a*s + 1);
W{6}  = K/(T2^2*s^2 + T1*s + 1);
W{7}  = K/(T2^2*s^2 + 1);
W{8}  = K/(s*(T*s + 1));
W{9}  = K*s/(T*s + 1);
W{10} = K*(T*s + 1);
W{11} = K*(T*s + 1)/s;

s_sym = sym('s');
t_sym = sym('t');

W_sym = cell(11,1);
W_sym{1}  = K/s_sym;
W_sym{2}  = K*s_sym;
W_sym{3}  = K;
W_sym{4}  = K/(T*s_sym + 1);
W_sym{5}  = K/(T2^2*s_sym^2 + T1a*s_sym + 1);
W_sym{6}  = K/(T2^2*s_sym^2 + T1*s_sym + 1);
W_sym{7}  = K/(T2^2*s_sym^2 + 1);
W_sym{8}  = K/(s_sym*(T*s_sym + 1));
W_sym{9}  = K*s_sym/(T*s_sym + 1);
W_sym{10} = K*(T*s_sym + 1);
W_sym{11} = K*(T*s_sym + 1)/s_sym;

fig1 = figure('Name', 'h(t)', 'Position', [50 50 1200 800]);
for i = 1:numel(names)
    subplot(3, 4, i); hold on; grid on;
    
    if i == 2
        y = zeros(size(time)); 
        y(1) = K / dt; 
        plot(time, y, 'b', 'LineWidth', 1);
    elseif i == 3
        y = K * ones(size(time));
        plot(time, y, 'b', 'LineWidth', 1);
    elseif i == 10
        y = K * ones(size(time));
        y(1) = K + (K * T) / dt;
        plot(time, y, 'b', 'LineWidth', 1);
    else
        [y, t_out] = step(W{i}, time);
        plot(t_out, y, 'b', 'LineWidth', 1);
    end
    
    Y_s = W_sym{i} * (1/s_sym);
    y_t = ilaplace(Y_s, s_sym, t_sym);
    y_analyt = double(subs(y_t, t_sym, time));
    plot(time, y_analyt, 'r--', 'LineWidth', 2);
    
    xlim([0 5])
    title(names{i}, 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('t, c'); ylabel('h(t)');
end

fig2 = figure('Name', 'w(t)', 'Position', [100 100 1400 800]);
for i = 1:numel(names)
    subplot(3, 4, i); hold on; grid on;
    
    if i == 2
        y = zeros(size(time));
        y(1) = K / (dt^2); 
        y(2) = -K / (dt^2);
        plot(time, y, 'r', 'LineWidth', 1);
    elseif i == 3
        y = zeros(size(time)); 
        y(1) = K / dt;
        plot(time, y, 'r', 'LineWidth', 1);
    elseif i == 10
        y = zeros(size(time)); 
        y(1) = (K * T) / (dt^2); 
        y(2) = -(K * T) / (dt^2);
        plot(time, y, 'r', 'LineWidth', 1);
    elseif i == 11
        y = K * ones(size(time)); 
        y(1) = K + (K * T) / dt;
        plot(time, y, 'r', 'LineWidth', 1);
        ylim([0 2*K]);
    else
        [y, t_out] = impulse(W{i}, time);
        plot(t_out, y, 'r', 'LineWidth', 1);
    end
    
    Y_s = W_sym{i};
    y_t = ilaplace(Y_s, s_sym, t_sym);
    y_analyt = double(subs(y_t, t_sym, time));
    plot(time, y_analyt, 'b--', 'LineWidth', 1.5);
    
    xlim([0 5])
    title(names{i}, 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('t, c'); ylabel('w(t)');
end

fig_var_KT = figure('Name', 'Варьирование K и T', 'Position', [100 100 1400 600]);
subplot(1, 2, 1);
hold on; grid on;
K_vals = [0.5, 1.0, 2.2, 3.5, 5.0];
for k_val = K_vals
    [y, t_out] = step(k_val / (T2^2 * s^2 + T1 * s + 1), time);
    plot(t_out, y, 'LineWidth', 1.5); 
end
title('Колебательное звено: варьирование K', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('t, c', 'FontSize', 10); 
ylabel('h(t)', 'FontSize', 10);
legend('K=0.5', 'K=1.0', 'K=2.2 (вар)', 'K=3.5', 'K=5.0', 'Location', 'northeast');

subplot(1, 2, 2);
hold on; grid on;
T_factors = [0.2, 0.5, 1.0, 2.0, 5.0];
for f = T_factors
    [y, t_out] = step(K / ((T * f) * s + 1), time);
    plot(t_out, y, 'LineWidth', 1.5); 
end
title('Апериодическое звено: варьирование T', 'FontSize', 11, 'FontWeight', 'bold');
xlabel('t, c', 'FontSize', 10); 
ylabel('h(t)', 'FontSize', 10);
legend('T/5 (0.04 c)', 'T/2 (0.1 c)', 'T (0.2 c)', '2T (0.4 c)', '5T (1.0 c)', 'Location', 'southeast');

fig5 = figure('Name', 'Влияние T1', 'Position', [100 100 900 500]);
hold on; grid on;
T1_vals = [0.5, 0.4, 0.25, 0.1, 0.02];
for t1_val = T1_vals
    [y, t_out] = step(K / (T2^2 * s^2 + t1_val * s + 1), time);
    plot(t_out, y, 'LineWidth', 1.5); 
end
title('Колебательное звено: уменьшение T1 (приближение к консервативному)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('t, c'); ylabel('h(t)');
legend('T1=0.5', 'T1=0.4 (вар)', 'T1=0.25', 'T1=0.1', 'T1=0.02', 'Location', 'northeast');
