clc;
clear; 
close all;
s = tf('s');
T_pp_req = 0.8;                
M_req = 1.35;
sigma_req = 30;
W_open = 200 / (s * (0.1*s + 1) * (0.02*s + 1) * (0.01*s + 1));

W_cl_open = feedback(W_open, 1);
poles_open_cl = pole(W_cl_open);
disp('Полюса исходной замкнутой системы:');
disp(poles_open_cl);

gamma_req = asind(1 / M_req);           
wc_des = ceil(2*pi / T_pp_req);             
Wc = evalfr(W_open, 1i * wc_des);
gamma_base = 180 + rad2deg(angle(Wc));
phi_max = gamma_req - gamma_base; 
alpha = (1 + sind(phi_max)) / (1 - sind(phi_max));
T2 = 1 / (sqrt(alpha) * wc_des);
T1 = alpha * T2;
L_base = 20*log10(abs(Wc));
L_lead = 10*log10(alpha);
k = 10^(-(L_base + L_lead) / 20);

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
