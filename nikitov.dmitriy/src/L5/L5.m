clc;
clear;
close all;

s = tf('s');
W = 900 / (s * (1 + 0.08*s) * (1 + 0.02*s));
H0 = feedback(W, 1);
[gm0, pm0, ~, wc0] = margin(W);
disp('Полюсы исходной замкнутой системы:');
disp(pole(H0));
fprintf('До коррекции: wc = %.3f рад/с, GM = %.2f дБ, PM = %.2f град\n', ...
    wc0, 20*log10(gm0), pm0);

%% Синтез опережающего звена
Mlimit = 1.5;
sigmaLimit = 37;
phaseTarget = asind(1/Mlimit);
b = 0.1 / tand(phaseTarget);
wcMax = (b + sqrt(b^2 + 4*0.0016)) / (2*0.0016);
bestTime = Inf;

for wc = 0.1:0.1:wcMax
    plantPhase = -90 - atand(0.08*wc) - atand(0.02*wc);
    lead = max(0, phaseTarget - (180 + plantPhase));
    alpha = (1 + sind(lead)) / (1 - sind(lead));
    Tz = sqrt(alpha) / wc;
    Tp = 1 / (wc * sqrt(alpha));
    K = 1 / (abs(evalfr(W, 1i*wc)) * sqrt(alpha));

    C = K * (1 + Tz*s) / (1 + Tp*s);
    H = feedback(C*W, 1);
    if ~isstable(H) || getPeakGain(H, 1e-7) > Mlimit
        continue;
    end
    info = stepinfo(H, 'SettlingTimeThreshold', 0.02);
    if info.Overshoot <= sigmaLimit && info.SettlingTime < bestTime
        bestTime = info.SettlingTime;
        best = [wc, K, Tz, Tp];
    end
end

wc = best(1);
K = best(2);
Tz = best(3);
Tp = best(4);
C = K * (1 + Tz*s) / (1 + Tp*s);
Wd = C * W;
H = feedback(Wd, 1);
[gm, pm] = margin(Wd);
[M, wPeak] = getPeakGain(H, 1e-7);
info = stepinfo(H, 'SettlingTimeThreshold', 0.02);
fprintf('wc = %.3f рад/с, K = %.9f, Tz = %.9f с, Tp = %.9f с\n', wc, K, Tz, Tp);
fprintf('Частоты сопряжения ПКУ: %.4f и %.4f рад/с\n', 1/Tz, 1/Tp);
fprintf('После коррекции: GM = %.2f дБ, PM = %.2f град, M = %.4f при %.3f рад/с\n', ...
    20*log10(gm), pm, M, wPeak);
fprintf('Tпп (2%%) = %.4f с, sigma = %.3f %%\n', info.SettlingTime, info.Overshoot);
disp('Полюсы скорректированной замкнутой системы:');
disp(pole(H));

%% Переходные процессы
figure('Name', 'Переходные процессы', 'Theme', 'light', 'Position', [70 100 920 475]);
tiledlayout(1, 2, 'Padding', 'compact');
nexttile;
t0 = 0:0.0002:0.23;
plot(t0, step(H0, t0), 'LineWidth', 1.5);
grid on;
title('До коррекции');
xlabel('t, с');
ylabel('y(t)');
nexttile;
t1 = 0:0.0002:0.30;
plot(t1, step(H, t1), 'LineWidth', 1.5);
yline(0.98, ':');
yline(1.02, ':');
grid on;
title('После коррекции');
xlabel('t, с');
ylabel('y(t)');

%% Исходные и желаемые характеристики
w = logspace(-1, 4, 1500);
[m0, p0] = bode(W, w);
[md, pd] = bode(Wd, w);
figure('Name', 'Разомкнутый контур', 'Theme', 'light', 'Position', [90 90 920 680]);
tiledlayout(2, 1, 'Padding', 'compact');
nexttile;
semilogx(w, 20*log10(m0(:)), w, 20*log10(md(:)), 'LineWidth', 1.5);
grid on;
yline(0, ':', 'HandleVisibility', 'off');
xline(wc, ':', 'HandleVisibility', 'off');
ylabel('L, дБ');
legend('Исходная', 'Желаемая', 'Location', 'southwest');
nexttile;
semilogx(w, p0(:), w, pd(:), 'LineWidth', 1.5);
grid on;
yline(-180, ':');
xline(wc, ':');
xlabel('\omega, рад/с');
ylabel('\phi, град');

%% Геометрическое вычитание асимптот
L0 = 20*log10(900./w) - 20*log10(max(1, 0.08*w)) - 20*log10(max(1, 0.02*w));
Ld = L0 + 20*log10(K) + 20*log10(max(1, Tz*w)) - 20*log10(max(1, Tp*w));
figure('Name', 'Асимптотические ЛАЧХ', 'Theme', 'light', 'Position', [110 100 920 475]);
semilogx(w, L0, '--', w, Ld, '-.', w, Ld-L0, 'LineWidth', 1.5);
grid on;
xline(1/Tz, ':', 'HandleVisibility', 'off');
xline(1/Tp, ':', 'HandleVisibility', 'off');
xlabel('\omega, рад/с');
ylabel('L, дБ');
legend('Исходная', 'Желаемая', 'ПКУ', 'Location', 'southwest');

%% Частотные характеристики замкнутой системы после коррекции
[m, phase] = bode(H, w);
figure('Name', 'Скорректированная замкнутая САУ', 'Theme', 'light', ...
    'Position', [130 90 920 680]);
tiledlayout(2, 1, 'Padding', 'compact');
nexttile;
semilogx(w, 20*log10(m(:)), 'LineWidth', 1.5);
hold on;
plot(wPeak, 20*log10(M), 'ko', 'MarkerFaceColor', 'k');
grid on;
xlim([1, 3000]);
ylabel('L, дБ');
nexttile;
semilogx(w, phase(:), 'LineWidth', 1.5);
grid on;
xlim([1, 3000]);
xlabel('\omega, рад/с');
ylabel('\phi, град');
