%% Аналитические вычисления
syms s t K;
% Передаточная функция динамического звена (в зависимости от варианта)
W = ???
% 1. Единичное ступенчатое воздействие
U_step = laplace(???);
% 2. Выходной сигнал
Y_step = ???;
% 3. Переходная функция
h = ilaplace(???);
% 4. Импульсное воздействие
U_impulse = laplace(???);
% 5. Выходной сигнал
Y_impulse = ???;
% 6. Импульсная переходная функция
w = ilaplace(Y_impulse, s, t);

%% 7. Построение переходных характеристик
K_values = 1:10;
time = 0:0.01:10;

figure;
hold on;
grid on;

for k = K_values
    h_k = ???;
    y = ???;
    plot(time, y, 'DisplayName', ['K = ', num2str(k)]);
end

title('Переходные характеристики интегрирующего звена');
xlabel('Время t, с');
ylabel('h(t)');
legend('show', 'Location', 'northwest');
hold off;

%% 8. Построение импульсных характеристик
figure;
hold on;
grid on;

for k = K_values
    w_k = ???;
    y = ???;
    plot(time, y, 'DisplayName', ['K = ', num2str(k)]);
end

title('Импульсные характеристики интегрирующего звена');
xlabel('Время t, с');
ylabel('w(t)');
legend('show', 'Location', 'northwest');
hold off;

%% 9. Проверка с помощью Control System Toolbox
s_tf = tf('s');

figure;
hold on;
grid on;

for k = K_values
    % Задаётся вручную
    W_check = ???;
    step(W_check, time);
end

title('Проверка переходной характеристики');
xlabel('Время t, с');
ylabel('Выходной сигнал');
hold off;

figure;
hold on;
grid on;
for k = K_values
    % Задаётся вручную
    W_check = ???;
    impulse(W_check, time);
end
title('Проверка импульсной характеристики');
xlabel('Время t, с');
ylabel('Выходной сигнал');
hold off;
