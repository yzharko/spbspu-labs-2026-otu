clc;
clear;
close all;

s = tf('s');
W0 = 300 / (s*(0.2*s + 1)*(0.01*s + 1));
G0 = feedback(W0, 1);

M_req = 1.3;
sigma_req = 27;
PM_target = 52.7;

wc_trials = 40:2.5:80;
best_settling = Inf;
wc = NaN;

for wc_trial = wc_trials
    [mag_trial, phase_trial] = bode(W0, wc_trial);
    mag_trial = squeeze(mag_trial);
    phase_trial = squeeze(phase_trial);
    phi_trial = mod(PM_target - (180 + phase_trial), 360);
    alpha_trial = (1 + sind(phi_trial))/(1 - sind(phi_trial));
    T2_trial = 1/(wc_trial*sqrt(alpha_trial));
    T1_trial = alpha_trial*T2_trial;
    Wlead_trial = (T1_trial*s + 1)/(T2_trial*s + 1);
    [mag_lead_trial, ~] = bode(Wlead_trial, wc_trial);
    mag_lead_trial = squeeze(mag_lead_trial);
    Kc_trial = 1/(mag_trial*mag_lead_trial);
    G_trial = feedback(W0*Kc_trial*Wlead_trial, 1);
    if ~isstable(G_trial)
        continue;
    end
    M_trial = getPeakGain(G_trial);
    S_trial = stepinfo(G_trial);
    if M_trial <= M_req && S_trial.Overshoot <= sigma_req && ...
            S_trial.SettlingTime < best_settling
        best_settling = S_trial.SettlingTime;
        wc = wc_trial;
    end
end

assert(~isnan(wc), 'Не найдена частота среза, удовлетворяющая требованиям');

[mag0, phase0] = bode(W0, wc);
mag0 = squeeze(mag0);
phase0 = squeeze(phase0);
phi_lead = mod(PM_target - (180 + phase0), 360);
alpha = (1 + sind(phi_lead))/(1 - sind(phi_lead));
T2 = 1/(wc*sqrt(alpha));
T1 = alpha*T2;
Wlead = (T1*s + 1)/(T2*s + 1);
[mag_lead, ~] = bode(Wlead, wc);
mag_lead = squeeze(mag_lead);
Kc = 1/(mag0*mag_lead);
Wku = Kc*Wlead;
W = W0*Wku;
G = feedback(W, 1);

[Gm, Pm, ~, ~] = margin(W);
Gm_dB = 20*log10(Gm);
M_fact = getPeakGain(G);
info = stepinfo(G);

fprintf('W_ку(s) = %.4f*(%.4f*s + 1)/(%.6f*s + 1)\n', Kc, T1, T2);
fprintf('Pm = %.2f°, Gm = %.2f дБ, M = %.4f, Tпп = %.4f с, sigma = %.1f%%\n', ...
    Pm, Gm_dB, M_fact, info.SettlingTime, info.Overshoot);

w = logspace(-1, 4, 3000);
asy = @(T) 20*log10(max(1, w*T));
L0 = 20*log10(300./w) - asy(0.2) - asy(0.01);
Ldes = 20*log10(300*Kc./w) - asy(0.2) - asy(0.01) + asy(T1) - asy(T2);
Lku = Ldes - L0;

[mag0_all, ph0_all] = bode(W0, w);
[magd_all, phd_all] = bode(W, w);

figure('Color', 'w');
subplot(2, 1, 1);
semilogx(w, 20*log10(squeeze(mag0_all)), 'b', ...
         w, 20*log10(squeeze(magd_all)), 'r', 'LineWidth', 1.5);
grid on;
ylabel('L, дБ');
legend('исходная', 'желаемая');
subplot(2, 1, 2);
semilogx(w, squeeze(ph0_all), 'b', ...
         w, squeeze(phd_all), 'r', 'LineWidth', 1.5);
grid on;
xlabel('\omega, рад/с');
ylabel('\psi, °');

figure('Color', 'w');
semilogx(w, L0, 'b', w, Ldes, 'r', w, Lku, 'g--', 'LineWidth', 1.5);
grid on;
xlabel('\omega, рад/с');
ylabel('L, дБ');
legend('L_0', 'L_ж', 'L_ку', 'Location', 'southwest');
title('Асимптотические ЛАЧХ');

figure('Color', 'w');
bode(G, w);
grid on;
title('ЛАЧХ и ЛФЧХ скорректированной замкнутой системы');

t0 = 0:0.0005:0.5;
t1 = 0:0.0005:0.2;

figure('Color', 'w');
[y0, t0] = step(G0, t0);
plot(t0, y0, 'LineWidth', 1.5);
grid on;
xlabel('t, с');
ylabel('y(t)');
title('Исходная система');

figure('Color', 'w');
[y1, t1] = step(G, t1);
plot(t1, y1, 'LineWidth', 1.5);
grid on;
xlabel('t, с');
ylabel('y(t)');
title('Скорректированная система');
