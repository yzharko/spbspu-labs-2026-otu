clc;
clear;
close all;

s = tf('s');

W = 100 / (s * (0.5*s + 1));

Wku = 0.4472 * (0.5*s + 1) / (0.025*s + 1);

Wcorr = Wku * W;

Wz = feedback(W, 1);
Wz_corr = feedback(Wcorr, 1);

figure;
step(Wz, Wz_corr);
grid on;
legend('Исходная система', 'Скорректированная система');
title('Переходные характеристики');

info_before = stepinfo(Wz);
info_after = stepinfo(Wz_corr);

fprintf('ИСХОДНАЯ СИСТЕМА\n');
fprintf('Tпп = %.4f с\n', info_before.SettlingTime);
fprintf('Перерегулирование = %.2f %%\n\n', info_before.Overshoot);

fprintf('СКОРРЕКТИРОВАННАЯ СИСТЕМА\n');
fprintf('Tпп = %.4f с\n', info_after.SettlingTime);
fprintf('Перерегулирование = %.2f %%\n\n', info_after.Overshoot);

figure;
bode(W, Wcorr);
grid on;
legend('Исходная', 'Желаемая (скорректированная)');
title('ЛАЧХ и ЛФЧХ исходной и желаемой систем');

figure;
bode(Wku);
grid on;
title('ЛАЧХ и ЛФЧХ корректирующего устройства');

figure;
margin(W);
grid on;
title('Исходная система');

figure;
margin(Wcorr);
grid on;
title('Скорректированная система');

[~, ~, ~, wc_before] = margin(W);
[~, ~, ~, wc_after] = margin(Wcorr);

fprintf('Частота среза исходной системы = %.4f рад/с\n', wc_before);
fprintf('Частота среза скорректированной системы = %.4f рад/с\n', wc_after);

wM = logspace(-2, 4, 10000);

mag_before = squeeze(abs(freqresp(Wz, wM)));
mag_after = squeeze(abs(freqresp(Wz_corr, wM)));

M_before = max(mag_before);
M_after = max(mag_after);

fprintf('\nПоказатель колебательности исходной системы M = %.4f\n', M_before);
fprintf('Показатель колебательности скорректированной системы M = %.4f\n', M_after);

w = logspace(-1, 3, 5000);

L0_asym = 20*log10(100 ./ w);
L0_asym = L0_asym - 20*max(0, log10(w / 2));

Lj_asym = 20*log10(44.72 ./ w);
Lj_asym = Lj_asym - 20*max(0, log10(w / 40));

Lku_asym = Lj_asym - L0_asym;

resp_ku = squeeze(freqresp(Wku, w));
Lku_exact = 20*log10(abs(resp_ku));

figure;

subplot(2,1,1);
semilogx(w, L0_asym, 'LineWidth', 1.5);
hold on;
semilogx(w, Lj_asym, 'LineWidth', 1.5);
yline(0, 'k:');
grid on;
legend('Исходная', 'Желаемая', 'Location', 'southwest');
ylabel('L(\omega), дБ');
title('Асимптотические ЛАЧХ исходной и желаемой систем');

subplot(2,1,2);
semilogx(w, Lku_asym, 'LineWidth', 1.5);
hold on;
semilogx(w, Lku_exact, '--', 'LineWidth', 1.5);
xline(2, ':');
xline(40, ':');
grid on;
legend('Геометрическая разность', 'Точная ЛАЧХ W_{ку}', 'Location', 'northwest');
xlabel('\omega, рад/с');
ylabel('L_{ку}(\omega), дБ');
title('Геометрическое вычитание ЛАЧХ');

resp_closed = squeeze(freqresp(Wz_corr, w));
mag_closed = abs(resp_closed);
phase_closed = unwrap(angle(resp_closed)) * 180 / pi;
L_closed = 20*log10(mag_closed);

[M_closed, idx] = max(mag_closed);
omega_res = w(idx);

figure;

subplot(2,1,1);
semilogx(w, L_closed, 'LineWidth', 1.5);
hold on;
plot(omega_res, 20*log10(M_closed), 'ro', 'MarkerFaceColor', 'r');
grid on;
ylabel('L(\omega), дБ');
title(sprintf('ЛАЧХ скорректированной замкнутой системы, M = %.4f', M_closed));

subplot(2,1,2);
semilogx(w, phase_closed, 'LineWidth', 1.5);
grid on;
xlabel('\omega, рад/с');
ylabel('\phi(\omega), град');
title('ЛФЧХ скорректированной замкнутой системы');

fprintf('Резонансная частота скорректированной замкнутой системы = %.4f рад/с\n', omega_res);
