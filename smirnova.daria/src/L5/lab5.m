clc;
clear; 
close all;
s = tf('s');
W_open = 200 / (s * (0.1*s + 1) * (0.02*s + 1) * (0.01*s + 1));
W_cl_open = feedback(W_open, 1); %замыкает систему 1 отриц. обратной связью

poles_open_cl = pole(W_cl_open);
disp('Полюса исходной замкнутой системы:');
disp(poles_open_cl);
k = 0.04;
T1 = 0.2;
T2 = 0.05;
W_cor = k * (T1*s + 1) / (T2*s + 1);
W_corr_open = W_cor * W_open;
W_cl_corr = feedback(W_corr_open, 1);
poles_corr_cl = pole(W_cl_corr);
info_corr = stepinfo(W_cl_corr);
disp('Полюса скорректированной замкнутой системы:');
disp(poles_corr_cl);
disp('Показатели переходного процесса скорректированной системы:');
disp(info_corr);

[mag_cl, ~, ~] = bode(W_cl_corr);%возвращает массивы амплитуды, фазы и частот замкнутой системы
M_fact = max(squeeze(mag_cl));
fprintf('Фактический показатель колебательности M = %.3f\n', M_fact);

figure('Name', 'Переходные характеристики');
subplot(2, 1, 1);
step(W_cl_open, 1.5);
grid on;
title('Исходная замкнутая система (неустойчива)');
xlabel('Время, с');
ylabel('Выходной сигнал');
subplot(2, 1, 2);
step(W_cl_corr, 1.5);
grid on;
title('Скорректированная замкнутая система');
xlabel('Время, с');
ylabel('Выходной сигнал');

figure('Name', 'ЛАЧХ и ЛФЧХ до и после коррекции');
bode(W_open, 'r--', W_corr_open, 'b-');
grid on;
legend('Исходная W_{open}(s)', 'Скорректированная W_{corr}(s)', 'Location', 'best');
title('ЛАЧХ и ЛФЧХ разомкнутых систем');

figure('Name', 'Характеристики корректирующего устройства');
bode(W_cor);
grid on;
title('ЛАЧХ и ЛФЧХ корректирующего устройства W_{cor}(s)');
