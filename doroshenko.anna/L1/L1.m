syms s t;
syms K ;
T = sym(0.6);  
T1 = sym(0.8);  
T2 = sym(0.7); 
T1_aperiodic = sym(1.5);

W_names = {"Интегрирующее", ...
    "Дифференцирующее", ...
    "Усилительное", ...
    "Апериодическое 1-го порядка", ...
    "Апериодическое 2-го порядка", ...
    "Колебательное", ...
    "Консервативное", ...
    "Инерционное интегрирующее", ...
    "Реальное дифференцирующее", ...
    "Форсирующее 1-го порядка", ...
    "Изодромное"};

W = {K/s, ... 
    K*s, ...                            
    K, ...                        
    K/(T*s+1), ...                      
    K/(T2^2*s^2 + T1_aperiodic*s + 1), ...      
    K/(T2^2*s^2 + T1*s + 1), ...        
    K/(T*s^2 + 1), ...                  
    K/(s*(T*s+1)), ...                  
    K*s/(T*s+1), ...                    
    K*(T*s+1), ...                      
    K*(T*s+1)/s};                       
    
U_step = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);

h = cell(size(W));
w = cell(size(W));

for i = 1:numel(W)
    Y_step = W{i} * U_step;
    h{i} = ilaplace(Y_step, s, t);
    Y_impulse = W{i} * U_impulse;
    w{i} = ilaplace(Y_impulse, s, t);

    fprintf('\n%s звено\n', W_names{i});
    disp(['h(t) = ', char(h{i})]);
    disp(['w(t) = ', char(w{i})]);
end

K_values = [0.9, 1:10];
time = 0:0.01:10;

for i = 1:numel(W)
    figure;
    hold on;
    grid on;
    for k = K_values
        h_k = subs(h{i}, K, k);
        y = double(subs(h_k, t, time));
        plot(time, y, 'DisplayName', ['K = ', num2str(k)]);
    end
    title(['Переходные характеристики:', W_names{i}]);
    xlabel('Время t, с');
    ylabel('h(t)');
    legend('show', 'Location', 'northwest');
    hold off;
end

for i = 1:numel(W)
    figure;
    hold on;
    grid on;
    for k = K_values
        w_k = subs(w{i}, K, k);
        y = double(subs(w_k, t, time));
        plot(time, y, 'DisplayName', ['K = ', num2str(k)]);
    end
    title(['Импульсные характеристики:', W_names{i}]);
    xlabel('Время t, с');
    ylabel('w(t)');
    legend('show', 'Location', 'northwest');
    hold off;
end

s_tf = tf('s');
K_ch = 0.9;  
T_ch = double(T);  
T1_ch = double(T1);
T2_ch = double(T2);  
T1a_ch = double(T1_aperiodic);

for i = 1:numel(W_names)
    W_ch = W_by_index(i, 1, T_ch, T1_ch, T2_ch, T1a_ch, s_tf);
    if isnumeric(W_ch)
        W_ch = tf(W_ch);
    end
    if ~isproper(W_ch)
        continue;
    end
    figure;
    hold on;
    grid on;
    for k=K_values
        W_ch = W_by_index(i, k, T_ch, T1_ch, T2_ch, T1a_ch, s_tf);
        if isnumeric(W_ch)
            W_ch = tf(W_ch);
        end
        if ~isproper(W_ch)
            continue;
        end
        [y, tt] = step(W_ch, time);
        plot(tt, y, 'DisplayName', ['K = ', num2str(k)]);
    end
    title(['Проверка переходной характеристики: ', W_names{i}]);
    xlabel('Время t, с');
    ylabel('Выходной сигнал');
    legend('show', 'Location', 'northwest');
    hold off;
end

for i = 1:numel(W_names)
    W_ch = W_by_index(i, 1, T_ch, T1_ch, T2_ch, T1a_ch, s_tf);
    if isnumeric(W_ch)
        W_ch = tf(W_ch);
    end
    if ~isproper(W_ch)
        continue;
    end
    figure;
    hold on;
    grid on;
    for k=K_values
        W_ch = W_by_index(i, k, T_ch, T1_ch, T2_ch, T1a_ch, s_tf);
        if isnumeric(W_ch)
            W_ch = tf(W_ch);
        end
        if ~isproper(W_ch)
            continue;
        end
        [y, tt] = impulse(W_ch, time);
        plot(tt, y, 'DisplayName', ['K = ', num2str(k)]);
    end
    title(['Проверка импульсной характеристики: ', W_names{i}]);
    xlabel('Время t, с');
    ylabel('Выходной сигнал');
    legend('show', 'Location', 'northwest');
    hold off;
end

T_list = T_ch * [1/5, 1/2, 1, 2, 5];

figure; 
hold on; 
grid on;
fprintf('\nАпериодическое звено, K = %.1f\n', K_ch);
fprintf('   T      время переходного процесса\n');
for j = 1:numel(T_list)
    Wj = K_ch / (T_list(j)*s_tf + 1);
    [y, tt] = step(Wj, 0:0.01:20);
    plot(tt, y, 'DisplayName', ['T = ', num2str(T_list(j))]);
    info = stepinfo(Wj);
    fprintf('%6.2f   %8.2f c\n', T_list(j), info.SettlingTime);
end
title('Апериодическое звено');
xlabel('Время t, с');
ylabel('h(t)');
legend('show', 'Location', 'southeast');
hold off;

T1_list = [0.8, 0.6, 0.4, 0.2, 0.1];

figure; 
hold on; 
grid on;
fprintf('\nКолебательное звено, T2 = %.1f\n', T2_ch);
fprintf('   T1    перерегулирование, %%    время затухания, c\n');
for j = 1:numel(T1_list)
    Wj = K_ch / (T2_ch^2*s_tf^2 + T1_list(j)*s_tf + 1);
    [y, tt] = step(Wj, 0:0.01:40);
    plot(tt, y, 'DisplayName', ['T1 = ', num2str(T1_list(j))]);
    info = stepinfo(Wj);
    fprintf('%5.2f    %10.2f              %10.2f\n', ...
        T1_list(j), info.Overshoot, info.SettlingTime);
end
title('Колебательное звено');
xlabel('Время t, с');
ylabel('h(t)');
legend('show', 'Location', 'northeast');
hold off;

function Wi = W_by_index(i, Kx, Tx, T1x, T2x, T1ax, sx)
    switch i
        case 1,  Wi = Kx / sx;
        case 2,  Wi = Kx * sx;
        case 3,  Wi = Kx;
        case 4,  Wi = Kx / (Tx*sx + 1);
        case 5,  Wi = Kx / (T2x^2*sx^2 + T1ax*sx + 1);
        case 6,  Wi = Kx / (T2x^2*sx^2 + T1x*sx + 1);
        case 7,  Wi = Kx / (Tx*sx^2 + 1);
        case 8,  Wi = Kx / (sx*(Tx*sx + 1));
        case 9,  Wi = Kx*sx / (Tx*sx + 1);
        case 10, Wi = Kx * (Tx*sx + 1);
        case 11, Wi = Kx * (Tx*sx + 1) / sx;
    end
end
