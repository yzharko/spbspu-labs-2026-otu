clc;
clear;
close all;

%% Исходная система при единичной отрицательной обратной связи
s = tf('s');
K = 100;
omega0 = [10 50 100];                
M_target = 1.3;
sigma_limit = 25;                   
T_limit = 0.8;                     
band = 0.02;                         

W0 = K/(s*(1+s/omega0(1))*(1+s/omega0(2))*(1+s/omega0(3)));
Phi0 = feedback(W0, 1);
poles0 = pole(Phi0);
stable0 = isstable(Phi0);

t = linspace(0, T_limit, 80001).';
h0 = step(Phi0, t);
figure('Color','w', 'Name','Исходный переходный процесс');
plot(t, h0, 'Color',[0.65 0.23 0.16], 'LineWidth',1.6);
hold on;
yline(1, 'k:', 'Уровень задания', 'HandleVisibility','off');
grid on;
xlabel('t, с'); ylabel('h_0(t)');
title('Отклик исходной замкнутой системы');

%% Желаемая разомкнутая система
zeta = sqrt((1-sqrt(1-1/M_target^2))/2);
b_initial = 4*K*zeta^2;              
b = fzero(@(value) resonantPeak(K, value)-M_target, [100 200]);
omega_d = b*[1 5 10];
[M, omega_r] = resonantPeak(K, b);
Wd = K/(s*(1+s/omega_d(1))*(1+s/omega_d(2))*(1+s/omega_d(3)));
[gm0, pm0, w_phase0, wc0] = margin(W0);
[gm1, pm1, w_phase1, wc1] = margin(Wd);

%% Исходные и желаемые частотные характеристики в одном масштабе
w = logspace(0, 4, 5000);
[L0, phase0] = logResponse(W0, w);
[Ld, phase_d] = logResponse(Wd, w);
figure('Color','w', 'Name','Частотные характеристики разомкнутых контуров');
subplot(2,1,1);
semilogx(w, L0, 'r-', w, Ld, 'b-', 'LineWidth',1.4); hold on;
yline(0, 'k:', 'HandleVisibility','off');
semilogx(wc0, 0, 'ro', wc1, 0, 'bo', 'HandleVisibility','off');
grid on;
ylabel('L(\omega), дБ');
title('ЛАЧХ исходной и желаемой разомкнутых систем');
legend('Исходная','Желаемая','Location','southwest');
subplot(2,1,2);
semilogx(w, phase0, 'r-', w, phase_d, 'b-', 'LineWidth',1.4); hold on;
yline(-180, 'k:', 'HandleVisibility','off');
semilogx(wc0, pm0-180, 'ro', wc1, pm1-180, 'bo', 'HandleVisibility','off');
grid on;
xlabel('\omega, рад/с'); ylabel('\phi(\omega), град');
title('ЛФЧХ исходной и желаемой разомкнутых систем');

%% Геометрическое вычитание асимптотических ЛАЧХ
L0_asym = asymptoticMagnitude(w, K, omega0);
Ld_asym = asymptoticMagnitude(w, K, omega_d);
Lku_asym = Ld_asym-L0_asym;

Wku = (1+s/omega0(1))*(1+s/omega0(2))*(1+s/omega0(3)) / ...
      ((1+s/omega_d(1))*(1+s/omega_d(2))*(1+s/omega_d(3)));
[Lku, ~] = logResponse(Wku, w);

wc_asym = [sqrt(K*omega0(1)); K];
Lku_high_db = 20*log10(prod(omega_d)/prod(omega0));
frequency_summary = table([omega0; omega_d], [wc0; wc1], wc_asym, ...
    [w_phase0; w_phase1], [pm0; pm1], 20*log10([gm0; gm1]), ...
    'VariableNames', {'Breaks_rad_s','Crossover_rad_s','AsymptoticCrossover_rad_s', ...
                      'PhaseCrossover_rad_s','PhaseMargin_deg','GainMargin_dB'}, ...
    'RowNames', {'Исходная','Желаемая'});

figure('Color','w', 'Name','Синтез ПКУ геометрическим вычитанием');
subplot(2,1,1);
semilogx(w, L0_asym, 'r-', w, Ld_asym, 'b-', 'LineWidth',1.4); hold on;
yline(0, 'k:', 'HandleVisibility','off');
grid on; ylabel('L(\omega), дБ');
title('Асимптотические ЛАЧХ разомкнутых систем');
legend('Исходная','Желаемая','Location','southwest');
subplot(2,1,2);
semilogx(w, Lku_asym, 'b-', w, Lku, 'k--', 'LineWidth',1.4); hold on;
for frequency = [omega0 omega_d]
    xline(frequency, ':', 'Color',[0.55 0.55 0.55], 'HandleVisibility','off');
end
grid on;
xlabel('\omega, рад/с'); ylabel('L_{ку}(\omega), дБ');
title('ЛАЧХ последовательного корректирующего устройства');
legend('Геометрическая разность','Точная ЛАЧХ ПКУ','Location','northwest');

%% Проверка скорректированной замкнутой систем
Phi_full = feedback(ss(Wku)*ss(W0), 1);
Phi1 = feedback(minreal(Wku*W0), 1);
poles_full = pole(Phi_full);
poles1 = pole(Phi1);
stable1 = isstable(Phi_full) && isstable(Phi1);

h1 = step(Phi1, t);
h_inf = dcgain(Phi1);
S = stepinfo(h1, t, h_inf, 'SettlingTimeThreshold',band);
[L_closed, phase_closed] = logResponse(Phi1, w);
[M_closed, omega_r_closed] = getPeakGain(Phi1, 1e-8);

quality_summary = table([stable0; stable1; true], [NaN; S.SettlingTime; T_limit], ...
    [NaN; S.Overshoot; sigma_limit], [NaN; M_closed; M_target], ...
    'VariableNames', {'Stable','SettlingTime_s','Overshoot_pct','M'}, ...
    'RowNames', {'До коррекции','После коррекции','Требование'});
disp(quality_summary);

figure('Color','w', 'Name','Скорректированный переходный процесс');
plot(t, h1, 'b-', 'LineWidth',1.6); hold on;
yline(h_inf*(1-band), 'k--'); yline(h_inf*(1+band), 'k--');
xline(S.SettlingTime, 'r:', sprintf('Tпп = %.4f с', S.SettlingTime));
plot(S.PeakTime, S.Peak, 'ro', 'MarkerFaceColor','r');
xlim([0 min(t(end), 2*S.SettlingTime)]); grid on;
xlabel('t, с'); ylabel('h_1(t)');
title('Отклик скорректированной замкнутой системы');

figure('Color','w', 'Name','Частотные характеристики после коррекции');
subplot(2,1,1);
semilogx(w, L_closed, 'b-', 'LineWidth',1.4); hold on;
semilogx(omega_r_closed, 20*log10(M_closed), 'ro', 'MarkerFaceColor','r');
grid on; ylabel('L(\omega), дБ');
title(sprintf('ЛАЧХ скорректированной замкнутой системы, M = %.2f', M_closed));
subplot(2,1,2);
semilogx(w, phase_closed, 'b-', 'LineWidth',1.4);
grid on;
xlabel('\omega, рад/с'); ylabel('\phi(\omega), град');
title('ЛФЧХ скорректированной замкнутой системы');

function [M, omega_r] = resonantPeak(K, b)
    rho = K/b;
    F = [0.0004, 0.0504, 1.05+0.04*rho, 1-2.6*rho, rho^2];
    x = fzero(@(x) polyval(polyder(F), x), [0 10]);
    M = rho/sqrt(polyval(F, x));
    omega_r = b*sqrt(x);
end

function [L, phase] = logResponse(sys, w)
    response = reshape(freqresp(sys, w), 1, []);
    L = 20*log10(abs(response));
    phase = unwrap(angle(response))*180/pi;
end

function L = asymptoticMagnitude(w, K, breaks)
    L = 20*log10(K./w);
    for frequency = breaks
        L = L-20*max(0, log10(w/frequency));
    end
end
