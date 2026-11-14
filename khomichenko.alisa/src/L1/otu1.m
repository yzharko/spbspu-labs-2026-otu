%% Аналитические вычисления
syms s t K T T1 T2 T1_5;
% Передаточная функция динамического звена (в зависимости от варианта)
W_list = {K/s, K*s, K, K/(T*s+1), K/(T2^2*s^2+T1_5*s+1), K/(T2^2*s^2+T1*s+1),...
    K/(T*s^2+1), K/(s*(T*s+1)), (K*s)/(T*s+1), K*(T*s+1), (K*(T*s+1))/s};

W_name = {"Интегрирующего ", "Дифференцирующего ", "Усилительного ",...
    "Апериодического 1-го порядка ", "Апериодического 2-го порядка ", ...
    "Колебательного ", "Консервативного ", "Инерционного интегрирующего ",...
    "Реального дифференцирующего ", "Форсирующего 1-го порядка ", "Изодромного "};

T_val = 0.4;
T1_val = 0.5;
T2_val = 0.7;
T1_5_val = 1.5;
s_tf = tf('s');

for i = 1:numel(W_list)
    W = W_list{i};
    % 1. Единичное ступенчатое воздействие
    U_step = laplace(heaviside(t), t, s);
    % 2. Выходной сигнал
    Y_step = W * U_step;
    % 3. Переходная функция
    h = ilaplace(Y_step, s, t);
    % 4. Импульсное воздействие
    U_impulse = laplace(dirac(t), t, s);
    % 5. Выходной сигнал
    Y_impulse = W * U_impulse;
    % 6. Импульсная переходная функция
    w = ilaplace(Y_impulse, s, t);
    
    %% 7. Построение переходных характеристик
    K_values = [0.7, 1:10];
    time = 0:0.01:10;
    
    figure;
    hold on;
    grid on;
    
    for k = K_values
        h_k = subs(h, [K T T1 T2 T1_5], [k T_val T1_val T2_val T1_5_val]);
        y = double(subs(h_k, t, time));
        plot(time, y, 'DisplayName', ['K = ', num2str(k)]);
    end

    if has(h, 'dirac')
        h_add_comm = ". Содержит δ-компоненту.";
    else
        h_add_comm = "";
    end

    title('Переходные характеристики ' + W_name{i} + "звена " + h_add_comm);
    xlabel('Время t, с');
    ylabel('h(t)');
    legend('show', 'Location', 'northwest');
    hold off;
    
    %% 8. Построение импульсных характеристик
    figure;
    hold on;
    grid on;
    
    for k = K_values
        w_k = subs(w, [K T T1 T2 T1_5], [k T_val T1_val T2_val T1_5_val]);
        y = double(subs(w_k, t, time));
        plot(time, y, 'DisplayName', ['K = ', num2str(k)]);
    end

    if has(w, 'dirac')
        w_add_comm = ". Содержит δ-компоненту.";
    else
        w_add_comm = "";
    end

    title('Импульсные характеристики ' + W_name{i} + "звена" + w_add_comm);
    xlabel('Время t, с');
    ylabel('w(t)');
    legend('show', 'Location', 'northwest');
    hold off;
      
    %% 9. Проверка с помощью Control System Toolbox
    W_check = getWcheck(i, 1, s_tf, T_val, T1_val, T2_val, T1_5_val);
    if isproper(W_check)
        
        figure;
        hold on;
        grid on;
        
        for k = K_values
            W_check = getWcheck(i, k, s_tf, T_val, T1_val, T2_val, T1_5_val);
            step(W_check, time);
        end
        
        title('Проверка переходной характеристики ' + W_name{i} + "звена" + h_add_comm);
        xlabel('Время t, с');
        ylabel('Выходной сигнал');
        hold off;
        
        figure;
        hold on;
        grid on;
        for k = K_values
            W_check = getWcheck(i, k, s_tf, T_val, T1_val, T2_val, T1_5_val);
            impulse(W_check, time);
        end
        title('Проверка импульсной характеристики ' + W_name{i} + "звена" + w_add_comm);
        xlabel('Время t, с');
        ylabel('Выходной сигнал');
        hold off;
    end
    %% изменение скорости протекания процесса при увеличении и уменьшении постоянной времени T в 2 и в 5 раз
    if i == 4
        T_values = [T_val/5, T_val/2, T_val, T_val*2, T_val*5];
        K_val = 0.7;

        figure;
        hold on;
        grid on;

        for j = 1:numel(T_values)

            h_T = subs(h, [K T], [K_val T_values(j)]);
            y = double(subs(h_T, t, time));

            plot(time, y, 'DisplayName', ...
                ['T = ', num2str(T_values(j))]);

        end

        title('Влияние постоянной времени T на переходный процесс');
        xlabel('Время t, с');
        ylabel('h(t)');
        legend('show', 'Location', 'southeast');
        hold off;

    end

    %% зафиксировать постоянную времени 𝑇2 и исследовать изменение перерегулирования и времени затухания при уменьшении коэффициента 𝑇1
    if i == 6
        T1_values = [T1_val, T1_val/2, T1_val/5];

        figure;
        hold on;
        grid on;
        time = 0:0.01:30
        for j = 1:numel(T1_values)

            h_osc = subs(h, [K T1 T2], [K_val T1_values(j) T2_val]);
            y = double(subs(h_osc, t, time));

            plot(time, y, 'DisplayName', ['T1 = ', num2str(T1_values(j))]);
        end

        title('Влияние коэффициента T1 на колебательное звено');
        xlabel('Время t, с');
        ylabel('h(t)');
        legend('show', 'Location', 'southeast');
        hold off;

    end
end


function W_check = getWcheck(i, k, s_tf, T_val, T1_val, T2_val, T1_5_val)

    switch i
        case 1
            W_check = k / s_tf;
        case 2
            W_check = k * s_tf;
        case 3
            W_check = tf(k);
        case 4
            W_check = k / (T_val*s_tf+1);
        case 5
            W_check = k/(T2_val^2*s_tf^2+T1_5_val*s_tf+1);
        case 6
            W_check = k/(T2_val^2*s_tf^2+T1_val*s_tf+1);
        case 7
            W_check = k/(T_val*s_tf^2+1);
        case 8
            W_check = k/(s_tf*(T_val*s_tf+1));
        case 9
            W_check = (k*s_tf)/(T_val*s_tf+1);
        case 10
            W_check = k*(T_val*s_tf+1);
        case 11
            W_check = (k*(T_val*s_tf+1))/s_tf;
    end

end