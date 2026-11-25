clear;
clc;
close all;

%% Исследование исходной САУ
num = 100;
den = [0.00002 0.0032 0.13 1 0];
omega0 = [10 50 100];
W_open = tf(num, den);
W_closed = feedback(W_open, 1);

char_eq = den + [0 0 0 0 num];
a4 = char_eq(1);
a3 = char_eq(2);
a2 = char_eq(3);
a1 = char_eq(4);
a0 = char_eq(5);
D1 = a3;
D2 = a3*a2 - a4*a1;
D3 = a3*a2*a1 - a3^2*a0 - a4*a1^2;
D4 = a0*D3;
stable0 = a4 > 0 && all([D1 D2 D3 D4] > 0);
[gm0, pm0, w_phase0, wc0] = margin(W_open);
disp(sprintf('Определители Гурвица: D1=%.6g  D2=%.6g  D3=%.6g  D4=%.6g', ...
    D1, D2, D3, D4));

%% Начальная оценка и выбранная частота среза
Tpp_req = 0.8;
M_req = 1.3;
sigma_req = 25;
band = 0.02;
gamma_req = asind(1/M_req);
wc_start = pi/Tpp_req;
wc = 10;

%% Параметры опережающего звена
phi0_wc = -90 - sum(atand(wc./omega0));
gamma_base = 180 + phi0_wc;
delta_phi = gamma_req - gamma_base;

alpha = (1+sind(delta_phi))/(1-sind(delta_phi));
T2 = 1/(wc*sqrt(alpha));
T1 = alpha*T2;

L0_wc = 20*log10(abs(evalfr(W_open, 1i*wc)));
L_lead_wc = 20*log10(abs((1+1i*wc*T1)/(1+1i*wc*T2)));
L_K = -L0_wc - L_lead_wc;
K = 10^(L_K/20);
W_ku = K*tf([T1 1],[T2 1]);

%% Скорректированная САУ и геометрическое вычитание ЛАЧХ
W_final = W_ku*W_open;
W_closed_final = feedback(W_final, 1);
stable1 = isstable(W_closed_final);
info = stepinfo(W_closed_final, 'SettlingTimeThreshold', band);
Tpp_fact = info.SettlingTime;
sigma_fact = info.Overshoot;
[M_fact, omega_res] = getPeakGain(W_closed_final, 1e-8);

omega_z = 1/T1;
omega_p = 1/T2;
[gm1, pm1, w_phase1, wc1] = margin(W_final);
M_from_pm = 1/sind(pm1);
disp(sprintf('Оценка M = 1/sin(Pm) = %.2f (требуется <= %.2f)', ...
    M_from_pm, M_req));
disp(sprintf('Tпп = %.3f с (предел %.1f с), sigma = %.1f %% (предел %.0f %%)', ...
    Tpp_fact, Tpp_req, sigma_fact, sigma_req));

w = logspace(-1, 3, 4000);
[L0, phase0] = logResponse(W_open, w);
[Ld, phase_d] = logResponse(W_final, w);
[Lku, phase_ku] = logResponse(W_ku, w);

L0_asym = plantAsymptote(w, num, omega0);
Ld_asym = plantAsymptote(w, num*K, [omega0 omega_p]) ...
    + 20*max(0, log10(w/omega_z));
Lku_asym = Ld_asym - L0_asym;

%% Графики исходной и скорректированной САУ
t0 = linspace(0, Tpp_req, 5001).';
h0 = step(W_closed, t0);
figure('Color','w','Name','Исходная замкнутая САУ');
plot(t0,h0,'r-','LineWidth',1.4); hold on;
yline(1,'k:','Уровень задания');
grid on; xlabel('t, с'); ylabel('h_0(t)');
title('Переходная характеристика исходной замкнутой САУ');

figure('Color','w','Name','Частотные характеристики разомкнутых контуров');
subplot(2,1,1);
semilogx(w,L0,'r-',w,Ld,'b-',w,Lku,'k--','LineWidth',1.3); hold on;
yline(0,'k:','HandleVisibility','off');
semilogx(wc0,0,'ro',wc1,0,'bo','HandleVisibility','off');
grid on; ylabel('L(\omega), дБ');
title('ЛАЧХ разомкнутых контуров и ПКУ');
legend('Исходная','Желаемая (скорректированная)','ПКУ','Location','southwest');
subplot(2,1,2);
semilogx(w,phase0,'r-',w,phase_d,'b-',w,phase_ku,'k--','LineWidth',1.3); hold on;
yline(-180,'k:');
semilogx(wc1,pm1-180,'bo',wc,delta_phi,'ko');
grid on; xlabel('\omega, рад/с'); ylabel('\phi(\omega), град');
title('ЛФЧХ разомкнутых контуров и ПКУ');

figure('Color','w','Name','Геометрическое вычитание ЛАЧХ');
subplot(2,1,1);
semilogx(w,L0_asym,'r-',w,Ld_asym,'b-','LineWidth',1.3); hold on;
yline(0,'k:','HandleVisibility','off');
grid on; ylabel('L(\omega), дБ');
title('Исходная и желаемая асимптотические ЛАЧХ');
legend('Исходная','Желаемая','Location','southwest');
subplot(2,1,2);
semilogx(w,Lku_asym,'b-',w,Lku,'k--','LineWidth',1.3); hold on;
xline(omega_z,':','1/T_1','HandleVisibility','off');
xline(omega_p,':','1/T_2','HandleVisibility','off');
grid on; xlabel('\omega, рад/с'); ylabel('L_{ку}(\omega), дБ');
title('ЛАЧХ последовательного корректирующего устройства');
legend('Геометрическая разность','Точная ЛАЧХ ПКУ','Location','southeast');

figure('Color','w','Name','Переходные характеристики до и после коррекции');
subplot(2,1,1);
plot(t0,h0,'r-',t0,step(W_closed_final,t0),'b-','LineWidth',1.4);
grid on; xlabel('t, с'); ylabel('h(t)');
title('Сравнение откликов в одном масштабе');
legend('До коррекции','После коррекции','Location','best');
subplot(2,1,2);
t1 = linspace(0,max(1.25*Tpp_req,1.25*Tpp_fact),5001).';
h_inf = dcgain(W_closed_final);
plot(t1,step(W_closed_final,t1),'b-','LineWidth',1.4); hold on;
yline(h_inf*(1-band),'k--');
yline(h_inf*(1+band),'k--');
xline(Tpp_fact,'r:',sprintf('Tпп = %.4f с',Tpp_fact));
plot(info.PeakTime,info.Peak,'ro','MarkerFaceColor','r');
grid on; xlabel('t, с'); ylabel('h_1(t)');
title('Скорректированный отклик и полоса установления ±2 %');

[L_closed, phase_closed] = logResponse(W_closed_final, w);
figure('Color','w','Name','Частотные характеристики скорректированной замкнутой САУ');
subplot(2,1,1);
semilogx(w,L_closed,'b-','LineWidth',1.4); hold on;
semilogx(omega_res,20*log10(M_fact),'ro','MarkerFaceColor','r');
grid on; ylabel('L(\omega), дБ');
title(sprintf('ЛАЧХ скорректированной замкнутой САУ, M = %.3f',M_fact));
subplot(2,1,2);
semilogx(w,phase_closed,'b-','LineWidth',1.4);
grid on; xlabel('\omega, рад/с'); ylabel('\phi(\omega), град');
title('ЛФЧХ скорректированной замкнутой САУ');

function [L, phase] = logResponse(sys, w)
    response = reshape(freqresp(sys, w), 1, []);
    L = 20*log10(abs(response));
    phase = unwrap(angle(response))*180/pi;
end

function L = plantAsymptote(w, gain, breaks)
    L = 20*log10(gain./w);
    for frequency = breaks
        L = L - 20*max(0, log10(w/frequency));
    end
end
