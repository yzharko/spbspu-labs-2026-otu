clear;
clc;
close all;

num = 40;
den = [0.04 1.2 1 0];
W_open = tf(num, den);
W_closed = feedback(W_open, 1);

char_eq = den + [0 0 0 num];
a0 = char_eq(1);
a1 = char_eq(2);
a2 = char_eq(3);
a3 = char_eq(4);
D1 = a1;
D2 = a1*a2 - a0*a3;
D3 = a3*D2;
fprintf('Определители Гурвица: D1 = %.2f, D2 = %.2f, D3 = %.2f\n', D1, D2, D3);
disp('Полюсы исходной замкнутой САУ:');
disp(pole(W_closed));

[Am, Pm, ~, Wcp] = margin(W_open);
fprintf('Исходная САУ: запас по амплитуде %.3f дБ, по фазе %.3f град., частота среза %.3f рад/с\n', ...
    20*log10(Am), Pm, Wcp);

figure('Name', 'Исходная САУ');
step(W_closed, 0:0.002:20);
grid on;
title('Переходная характеристика исходной замкнутой САУ');

wc = 150;
gamma_target = 60;

wj = 1i*wc;
mag_i = abs(evalfr(W_open, wj));
ph_i = -90 - rad2deg(angle(polyval(den(1:3), wj)));
gamma_base = 180 + ph_i;
delta_phi = gamma_target - gamma_base;
fprintf('Фаза исходной САУ на частоте среза %.2f град., необходимое опережение %.2f град.\n', ...
    ph_i, delta_phi);

phi_stage = delta_phi/2;
alpha = (1 + sind(phi_stage)) / (1 - sind(phi_stage));
T2 = 1/(wc*sqrt(alpha));
T1 = alpha*T2;
w_zero = 1/T1;
w_pole = 1/T2;

L_wc = 20*log10(mag_i);
L_stage_wc = 20*log10(sqrt(1+(T1*wc)^2) / sqrt(1+(T2*wc)^2));
L_K = -L_wc - 2*L_stage_wc;

plant_breaks = sort(abs(roots(den(1:3))));
w = sort([logspace(-1, 4, 1500), plant_breaks.', w_zero, w_pole, wc]);

L_open_asym = 20*log10(num) - 20*log10(w) ...
    - 20*log10(max(1, w/plant_breaks(1))) ...
    - 20*log10(max(1, w/plant_breaks(2)));
L_des_asym = L_open_asym + L_K ...
    + 40*log10(max(1, w/w_zero)) - 40*log10(max(1, w/w_pole));
L_ku_asym = L_des_asym - L_open_asym;

K = 10^(L_ku_asym(1)/20);
W_ku = K * tf([T1 1], [T2 1])^2;
fprintf('Исходные сопрягающие частоты: %.3f и %.3f рад/с\n', plant_breaks(1), plant_breaks(2));
fprintf('ПКУ: K = %.6f, T1 = %.8f с, T2 = %.8f с\n', K, T1, T2);
fprintf('Частоты излома ПКУ: %.3f и %.3f рад/с\n', w_zero, w_pole);
fprintf('W_ku(s) = %.6f * (%.8f*s + 1)^2 / (%.8f*s + 1)^2\n', K, T1, T2);

W_final = W_ku * W_open;
W_closed_final = feedback(W_final, 1);
[Am2, Pm2, ~, Wcp2] = margin(W_final);

disp('Полюсы скорректированной замкнутой САУ:');
disp(pole(W_closed_final));
fprintf('После коррекции: запас по амплитуде %.3f дБ, по фазе %.3f град., частота среза %.3f рад/с\n', ...
    20*log10(Am2), Pm2, Wcp2);

info = stepinfo(W_closed_final, 'SettlingTimeThreshold', 0.02);
Tpp_fact = info.SettlingTime;
sigma_fact = info.Overshoot;
[M_fact, w_peak] = getPeakGain(W_closed_final, 1e-5);
fprintf('Пик АЧХ замкнутой САУ на частоте %.3f рад/с\n', w_peak);

fprintf('Время установления: %.4f с\n', Tpp_fact);
fprintf('Перерегулирование: %.2f %%\n', sigma_fact);
fprintf('Показатель колебательности M: %.4f\n', M_fact);

[mag_open, phase_open] = bode(W_open, w);
mag_open = reshape(mag_open, 1, []);
phase_open = reshape(phase_open, 1, []);
L_open = 20*log10(mag_open);
L_des = L_open + L_K + 40*log10(hypot(1, w/w_zero)./hypot(1, w/w_pole));
phase_des = phase_open + 2*(atand(w/w_zero) - atand(w/w_pole));

figure('Name', 'Исходная и желаемая ЛАЧХ и ЛФЧХ');
subplot(2,1,1);
semilogx(w, L_open, w, L_des, 'LineWidth', 1.2);
grid on;
ylabel('L, дБ');
legend('Исходная', 'Желаемая (после коррекции)', 'Location', 'best');
title('ЛАЧХ разомкнутой САУ');
subplot(2,1,2);
semilogx(w, phase_open, w, phase_des, 'LineWidth', 1.2);
grid on;
xlabel('omega, рад/с');
ylabel('Фаза, град.');
legend('Исходная', 'Желаемая (после коррекции)', 'Location', 'best');
title('ЛФЧХ разомкнутой САУ');

figure('Name', 'Геометрическое вычитание ЛАЧХ');
semilogx(w, L_open_asym, '--', w, L_des_asym, '--', w, L_ku_asym, 'LineWidth', 1.2);
grid on;
xlabel('omega, рад/с');
ylabel('L, дБ');
legend('Исходная L_{исх}', 'Желаемая L_ж', 'ПКУ: L_ж - L_{исх}', 'Location', 'best');
title('Асимптотические ЛАЧХ');

figure('Name', 'Скорректированная замкнутая САУ');
bode(W_closed_final, w);
grid on;
title('ЛАЧХ и ЛФЧХ скорректированной замкнутой САУ');

figure('Name', 'Сравнение переходных характеристик');
step(W_closed, W_closed_final, 0:0.0001:0.5);
grid on;
legend('Исходная САУ', 'Скорректированная САУ', 'Location', 'best');
title('Переходные характеристики замкнутых САУ');

figure('Name', 'Переходный процесс после коррекции');
[y_final, t_final] = step(W_closed_final, 0:0.0001:0.25);
plot(t_final, y_final, 'LineWidth', 1.2);
grid on;
xlabel('Время, с');
ylabel('h(t)');
hold on;
yline(0.98, '--');
yline(1.02, '--');
title('Скорректированная САУ и полоса ±2%');
