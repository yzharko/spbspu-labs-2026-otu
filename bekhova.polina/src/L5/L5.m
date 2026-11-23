clear;
clc;
close all;

s = tf('s');

num = 900;
den = [0.0016 0.1 1 0];

W_open = tf(num, den);
W_closed = feedback(W_open, 1);

char_eq = den + [0 0 0 num];

a3 = char_eq(1);
a2 = char_eq(2);
a1 = char_eq(3);
a0 = char_eq(4);

D1 = a2;
D2 = a2*a1 - a3*a0;
D3 = a0*D2;

fprintf('Определители Гурвица: D1=%.4f  D2=%.4f  D3=%.4f\n', D1, D2, D3);

poles_initial = pole(W_closed);

if all(real(poles_initial) < 0)
    fprintf('Исходная САУ устойчива\n');
else
    fprintf('Исходная САУ неустойчива\n');
end

[Am, Pm, Wcg, Wcp] = margin(W_open);

figure('Name','Исходная САУ');
subplot(2,1,1);
step(W_closed);
title('Переходная функция исходной САУ');
grid on;

subplot(2,1,2);
bode(W_open);
grid on;
title('ЛАЧХ/ЛФЧХ исходной разомкнутой САУ');

M_req = 1.5;
sigma_req = 37;

gamma_req = asind(1/M_req);

wc_values = linspace(10, 77, 250);

best_wc = NaN;
best_W_ku = [];
best_W_final = [];
best_W_closed_final = [];
best_T1 = NaN;
best_T2 = NaN;
best_K = NaN;
best_M = NaN;
best_sigma = NaN;
best_Pm = NaN;

for i = 1:length(wc_values)

    wc = wc_values(i);

    [mag_i, ph_i] = bode(W_open, wc);
    mag_i = squeeze(mag_i);
    ph_i = squeeze(ph_i);

    gamma_base = 180 + ph_i;
    delta_phi = gamma_req - gamma_base;

    if delta_phi <= 0 || delta_phi >= 89.5
        continue;
    end

    alpha = (1 + sind(delta_phi)) / (1 - sind(delta_phi));

    T2 = 1 / (wc * sqrt(alpha));
    T1 = alpha * T2;

    L_wc = 20*log10(mag_i);
    L_zv_wc = 20*log10(sqrt(1+(T1*wc)^2) / sqrt(1+(T2*wc)^2));

    L_K = -L_wc - L_zv_wc;
    K = 10^(L_K/20);

    W_ku = K * tf([T1 1], [T2 1]);

    W_final = W_ku * W_open;
    W_closed_final = feedback(W_final, 1);

    poles_final = pole(W_closed_final);

    if ~all(real(poles_final) < 0)
        continue;
    end

    [Am2, Pm2, Wcg2, Wcp2] = margin(W_final);

    if isnan(Pm2) || Pm2 <= 0 || Pm2 >= 90
        continue;
    end

    [M_fact, W_M] = getPeakGain(W_closed_final);

    info = stepinfo(W_closed_final);
    sigma_fact = info.Overshoot;

    if M_fact <= M_req && sigma_fact <= sigma_req
        best_wc = wc;
        best_W_ku = W_ku;
        best_W_final = W_final;
        best_W_closed_final = W_closed_final;
        best_T1 = T1;
        best_T2 = T2;
        best_K = K;
        best_M = M_fact;
        best_sigma = sigma_fact;
        best_Pm = Pm2;
    end
end

if isnan(best_wc)
    error('Не удалось подобрать ПКУ. Необходимо изменить диапазон частот.');
end

W_ku = best_W_ku;
W_final = best_W_final;
W_closed_final = best_W_closed_final;

T1 = best_T1;
T2 = best_T2;
K = best_K;

[M_fact, W_M] = getPeakGain(W_closed_final);

fprintf('\nПараметры желаемой системы:\n');
fprintf('Требуемый запас по фазе = %.2f град\n', gamma_req);
fprintf('Выбранная частота среза wc = %.4f рад/с\n', best_wc);

fprintf('\nПараметры корректирующего устройства:\n');
fprintf('T1 = %.6f\n', T1);
fprintf('T2 = %.6f\n', T2);
fprintf('K = %.6f\n', K);
fprintf('alpha = %.6f\n', T1/T2);

disp('Передаточная функция ПКУ:');
W_ku

disp('Скорректированная разомкнутая система:');
W_final

disp('Скорректированная замкнутая система:');
W_closed_final

[Am2, Pm2, Wcg2, Wcp2] = margin(W_final);

info = stepinfo(W_closed_final);

Tpp_fact = info.SettlingTime;
sigma_fact = info.Overshoot;

fprintf('\nПоказатели скорректированной системы:\n');
fprintf('Запас по амплитуде = %.4f\n', Am2);
fprintf('Запас по фазе = %.4f град\n', Pm2);
fprintf('Фактическое M = %.4f\n', M_fact);
fprintf('Частота максимума M = %.4f рад/с\n', W_M);
fprintf('Требуемое M = %.2f\n', M_req);
fprintf('Фактическое время переходного процесса = %.4f с\n', Tpp_fact);
fprintf('Фактическое перерегулирование = %.2f %%\n', sigma_fact);

figure('Name','ЛАЧХ и ЛФЧХ систем');
bode(W_open, W_ku, W_final);
grid on;
legend('Исходная W','ПКУ W_{ку}','Скорректированная W_{скор}','Location','best');
title('ЛАЧХ и ЛФЧХ исходной, корректирующей и скорректированной систем');

figure('Name','Переходные характеристики');
step(W_closed, W_closed_final);
grid on;
legend('Исходная САУ','Скорректированная САУ','Location','best');
title('Сравнение переходных характеристик');
