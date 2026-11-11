clear; 
clc; 
close all;
syms s t K_sym T_sym;
K  = 2.4;
T  = 0.4;
T1 = 0.8;
T2 = 0.9;
time = 0:0.01:10;
T1_ap = 2 * T2; %мин значение для пункта 5
eps_val = 0.001;
W{1}  = tf(K, [1 0]);                names{1}  = '1. Интегрирующее';
W{2} = tf([K 0], [eps_val 1]);       names{2}  = '2. Дифференцирующее (непричинное)';
W{3}  = tf(K, 1);                    names{3}  = '3. Усилительное';
W{4}  = tf(K, [T 1]);                names{4}  = '4. Апериодическое 1-го пор.';
W{5} = tf(K, [T2^2, T1_ap, 1]);      names{5}  = '5. Апериодическое 2-го пор.';
W{6}  = tf(K, [T2^2, T1, 1]);        names{6}  = '6. Колебательное';
W{7}  = tf(K, [T, 0, 1]);            names{7}  = '7. Консервативное';
W{8}  = tf(K, [T, 1, 0]);            names{8}  = '8. Интегр. с запаздыванием';
W{9}  = tf([K 0], [T 1]);            names{9}  = '9. Реальное дифференцирующее';
W{10} = tf([K*T K], [eps_val 1]);    names{10} = '10. Форсирующее (непричинное)';
W{11} = tf([K*T K], [1 0]);          names{11} = '11. Изодромное';

figure('Name', 'Переходные характеристики h(t)', 'Position', [100, 100, 1300, 900]);
length_of_w = length(names);
for i = 1:length_of_w
    subplot(4, 3, i);
    step(W{i}, time);
    if i == 2 || i == 10
        ylim([0 5]);
    end
    title(names{i}, 'FontSize', 9, 'FontWeight', 'bold');
    xlabel('t, c'); ylabel('h(t)'); grid on;
end
sgtitle('Переходные характеристики всех звеньев', 'FontSize', 14, 'FontWeight', 'bold');

figure('Name', 'Импульсные характеристики w(t)', 'Position', [150, 150, 1300, 900]);
for i = 1:length_of_w
    subplot(4, 3, i);
    impulse(W{i}, time);
    if i == 2 || i == 10 %Ограничили, чтобы не сжимать график
        ylim([-5 5]); 
    end
    title(names{i}, 'FontSize', 9, 'FontWeight', 'bold');
    xlabel('t, c'); ylabel('w(t)'); grid on;
end
sgtitle('Импульсные характеристики всех звеньев', 'FontSize', 14, 'FontWeight', 'bold');

figure('Name', 'Задание * : Исследование T для звена №4');
hold on; grid on;
T_factors = [1/5, 1/2, 1, 2, 5];
for f = T_factors
    T_curr = T * f;
    W_temp = tf(K, [T_curr 1]);
    step(W_temp, time);
end
title('Апериодическое звено №4: изменение T в 2 и 5 раз');
xlabel('Время t, с'); ylabel('h(t)');
legend('T/5 (0.08c)', 'T/2 (0.2c)', 'T (0.4c)', '2T (0.8c)', '5T (2.0c)');
hold off;

figure('Name', 'Задание * : Исследование T1 для звена №6');
hold on; grid on;
T1_values = [0.8, 0.5, 0.2, 0.05]; %Уменьшаем T1 к нулю
for t1_val = T1_values
    W_temp = tf(K, [T2^2, t1_val, 1]);
    step(W_temp, time);
end
title('Колебательное звено №6: приближение T1 к 0 (к консервативному)');
xlabel('Время t, с'); ylabel('h(t)');
legend('T1 = 0.8', 'T1 = 0.5', 'T1 = 0.2', 'T1 = 0.05');
hold off;

fprintf('\nПРОВЕРКА РЕЗУЛЬТАТОВ\n');
%Вектор символьных передаточных функций
syms s_chk t_chk;
W_sym_list = [
    K / s_chk;
    (K * s_chk) / (eps_val * s_chk + 1);
    K;
    K / (T * s_chk + 1);
    K / (T2^2 * s_chk^2 + T1_ap * s_chk + 1);
    K / (T2^2 * s_chk^2 + T1 * s_chk + 1);
    K / (T * s_chk^2 + 1);
    K / (s_chk * (T * s_chk + 1));
    (K * s_chk) / (T * s_chk + 1);
    (K * T * s_chk + K) / (eps_val * s_chk + 1);
    (K * T * s_chk + K) / s_chk
];

for idx = 1:length_of_w
    %Аналитический вывод h(t)
    Y_sym = W_sym_list(idx) * (1 / s_chk); %Умножаем передаточную функцию на входной сигнал
    h_func_sym = ilaplace(Y_sym, s_chk, t_chk);
    h_num = matlabFunction(h_func_sym, 'Vars', t_chk);
    y_anal = h_num(time);
    if isscalar(y_anal), y_anal = repmat(y_anal, size(time)); end

    [y_cst, ~] = step(W{idx}, time); %численный метод

    err = max(abs(y_anal(:) - y_cst(:)));
    fprintf('Звено №%-2d (%-30s): погрешность = %.2e\n', idx, names{idx}, err);
end
