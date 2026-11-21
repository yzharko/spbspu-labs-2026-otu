clear; clc; close all;

num = 1500;
den = [0.001  0.07  1  0];
W_open = tf(num, den);
W_closed = feedback(W_open, 1);

char_eq = [0.001  0.07  1  1500];
a3 = char_eq(4);
a2 = char_eq(3);
a1 = char_eq(2);
a0 = char_eq(1);

D1 = a2;
D2 = a2*a1 - a3*a0;
D3 = a3 * D2;

fprintf('Определители Гурвица: D1=%.2f  D2=%.2f  D3=%.2f\n', D1, D2, D3);

[Am, Pm, Wcg, Wcp] = margin(W_open);
fprintf('Запас по фазе исходной: %.2f град\n', Pm);
fprintf('Частота среза исходной: %.2f рад/с\n', Wcp);

fig1 = figure('Color','w');
subplot(2,1,1); step(W_closed, 0.2); grid on;
title('Переходная функция исходной САУ');
subplot(2,1,2); bode(W_open); grid on;
title('ЛАЧХ/ЛФЧХ исходной разомкнутой САУ');
saveas(fig1, 'Fig1_Original.png');

M_req = 1.5;
sigma_req = 37;

gamma_req = asind(1/M_req);
fprintf('Требуемый запас по фазе: %.2f град\n', gamma_req);

wc = 62;

[mag_i, ph_i] = bode(W_open, wc);
mag_i = squeeze(mag_i);
ph_i = squeeze(ph_i);
gamma_base = 180 + ph_i;
fprintf('Базовая фаза на wc: %.2f град\n', gamma_base);

delta_phi = gamma_req - gamma_base;
fprintf('Требуемое опережение: %.2f град\n', delta_phi);

alpha = (1 + sind(delta_phi)) / (1 - sind(delta_phi));
fprintf('Коэффициент alpha: %.2f\n', alpha);

T2 = 1 / (wc * sqrt(alpha));
T1 = alpha * T2;
fprintf('T1 = %.4f c, T2 = %.4f c\n', T1, T2);

L_wc = 20*log10(mag_i);
L_zv_wc = 20*log10(sqrt(1+(T1*wc)^2) / sqrt(1+(T2*wc)^2));
L_K = 0 - L_wc - L_zv_wc;
K = 10^(L_K/20);
fprintf('Коэффициент K = %.4f\n', K);

W_ku = K * tf([T1 1], [T2 1]);
disp('Передаточная функция ПКУ:');
W_ku

W_final = W_ku * W_open;
W_closed_final = feedback(W_final, 1);

[Am2, Pm2, Wcg2, Wcp2] = margin(W_final);
fprintf('Запас по фазе скорр.: %.2f град\n', Pm2);
fprintf('Частота среза скорр.: %.2f рад/с\n', Wcp2);

wM = logspace(-1, 3, 5000);
mag_closed = squeeze(abs(freqresp(W_closed_final, wM)));
[M_fact, ~] = max(mag_closed);
fprintf('Фактическое M = %.3f\n', M_fact);

info = stepinfo(W_closed_final);
fprintf('Tпп = %.4f с\n', info.SettlingTime);
fprintf('Перерегулирование = %.2f %%\n', info.Overshoot);

yesNo = @(b) subsref({'НЕТ','ДА'}, struct('type','{}','subs',{{b+1}}));
fprintf('Соответствие: sigma %.1f <= %.0f ? %s\n', info.Overshoot);

fig2 = figure('Color','w');
bode(W_open, W_ku, W_final); grid on;
legend('Исходная','ПКУ','Скорректированная','Location','best');
title('ЛАЧХ/ЛФЧХ — исходная, ПКУ, скорректированная');
saveas(fig2, 'Fig2_LACH_All.png');

fig3 = figure('Color','w');
step(W_closed_final, 0.5); grid on;
title('Переходная характеристика СКОРРЕКТИРОВАННОЙ САУ');
saveas(fig3, 'Fig3_Step_Corrected.png');

fig4 = figure('Color','w');
bode(W_closed_final); grid on;
title('ЛАЧХ/ЛФЧХ скорректированной замкнутой системы');
saveas(fig4, 'Fig4_Bode_Closed.png');

