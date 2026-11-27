clear;
clc;
close all;

%% 1. Исследование исходной САУ
num = 282;
den = [0.1 1 0];
W_open = tf(num, den);
W_closed = feedback(W_open, 1);
fprintf('Исходная САУ: определители Гурвица %.4g, %.4g\n', den(2), den(2)*num);

t = 0:0.0001:0.9;
y0 = squeeze(step(W_closed, t));
newFigure('Исходная САУ', 330);
plot(t, y0, 'b'); hold on;
yline(1, ':k'); yline([0.98 1.02], '--k');
grid on; xlabel('t, с'); ylabel('h(t)');

%% 2. Расчёт параметров желаемой ЛАЧХ
M_req = 1.5;
sigma_req = 37;
gamma_req = asind(1/M_req);
t1 = 0:0.000005:0.12;

fprintf('Подбор ωс при запасе по фазе %.6f°:\n', gamma_req);
for w_trial = [250 220 200 206 205]
    W_trial = feedback(lead(W_open, w_trial, gamma_req)*W_open, 1);
    trial_info = stepinfo(squeeze(step(W_trial, t1)), t1, 1, ...
        'SettlingTimeThreshold', 0.02);
    fprintf('ωс = %g: M = %.8f, σ = %.5f%%, Tпп = %.8f с\n', w_trial, ...
        getPeakGain(W_trial, 1e-9), trial_info.Overshoot, trial_info.SettlingTime);
end
wc = 205;

%% 3. Геометрическое вычитание и расчёт ПКУ
[W_ku, T1, T2, K, alpha, delta_phi, L_base, L_lead] = lead(W_open, wc, gamma_req);
Tgeom = sqrt(T1*T2);
assert(abs(1/Tgeom - wc) < 1e-10 && ...
    abs(atand((T1 - T2)/(2*Tgeom)) - delta_phi) < 1e-10);
w1 = 1/T1;
w2 = 1/T2;
w = sort([logspace(0, 4, 1000) 1/den(1) w1 wc w2]);
L0 = 20*log10(num./w) - 20*log10(max(1, w*den(1)));
Ld = 20*log10(K*num./w) - 20*log10(max(1, w*den(1))) ...
    + 20*log10(max(1, w/w1)) - 20*log10(max(1, w/w2));
Lku = Ld - L0;
fprintf('φmax = %.6f°, α = %.9f, Lисх(ωс) = %.6f дБ, Lоз(ωс) = %.6f дБ\n', ...
    delta_phi, alpha, L_base, L_lead);
fprintf('ПКУ: K = %.10f, T1 = %.11f с, T2 = %.11f с, ω1 = %.5f, ω2 = %.5f рад/с\n', ...
    K, T1, T2, w1, w2);

%% 4. Скорректированная САУ
W_final = W_ku*W_open;
W_closed_final = feedback(W_final, 1);
systems = {W_open, W_final; W_closed, W_closed_final};

%% 5. Показатели качества переходных процессов
time = {t, t1};
for i = 1:2
    [Gm, Pm, ~, wc_actual] = margin(systems{1,i});
    [M, wm] = getPeakGain(systems{2,i}, 1e-9);
    y = squeeze(step(systems{2,i}, time{i}));
    info = stepinfo(y, time{i}, 1, 'SettlingTimeThreshold', 0.02);
    fprintf('\nСАУ %d: ωс = %.6f рад/с, Δφ = %.6f°, ΔL = %.6f дБ\n', ...
        i, wc_actual, Pm, 20*log10(Gm));
    fprintf('M = %.8f при ω = %.6f рад/с, Tпп = %.8f с, σ = %.5f%%\n', ...
        M, wm, info.SettlingTime, info.Overshoot);
    disp(pole(systems{2,i}));
end
assert(isstable(W_closed_final) && M <= M_req && info.Overshoot <= sigma_req ...
    && abs(abs(evalfr(W_final, 1i*wc)) - 1) < 1e-10);

%% 6. Частотные и переходные характеристики
newFigure('Исходная и желаемая САУ', 500);
subplot(2,1,1);
semilogx(w, magnitude(W_open, w), 'b', w, magnitude(W_final, w), 'r'); hold on;
semilogx(w, L0, '--b', w, Ld, '--r');
yline(0, ':k'); grid on; ylabel('L, дБ');
legend('Исходная', 'Желаемая', 'Исходная асимпт.', 'Желаемая асимпт.', ...
    'Location', 'southwest');
subplot(2,1,2);
semilogx(w, phase(W_open, w), 'b', w, phase(W_final, w), 'r');
grid on; xlabel('\omega, рад/с'); ylabel('\phi, град');

newFigure('Геометрическое вычитание', 355);
semilogx(w, Lku, '-b', w, magnitude(W_ku, w), '--r'); hold on;
xline(w1, ':k', '\omega_1', 'LabelOrientation', 'horizontal');
xline(w2, ':k', '\omega_2', 'LabelOrientation', 'horizontal');
grid on; xlabel('\omega, рад/с'); ylabel('L_{ку}, дБ');
legend('Разность асимптот', 'Точная ЛАЧХ ПКУ', 'Location', 'southeast');

newFigure('Переходные характеристики', 465);
subplot(2,1,1);
plot(t, y0, 'b', t, squeeze(step(W_closed_final,t)), 'r'); hold on;
yline(1, ':k'); grid on; xlim([0 t(end)]);
ylabel('h(t)'); xlabel('t, с');
legend('Исходная', 'Скорректированная');
subplot(2,1,2);
plot(t1, squeeze(step(W_closed_final,t1)), 'r'); hold on;
yline(1, ':k'); yline([0.98 1.02], '--k');
xline(info.SettlingTime, ':k', 'T_{пп}', 'LabelOrientation', 'horizontal', ...
    'LabelVerticalAlignment', 'bottom');
grid on; xlim([0 0.08]); ylim([0 1.45]);
xlabel('t, с'); ylabel('h_{скор}(t)');

newFigure('Замкнутая скорректированная САУ', 420);
subplot(2,1,1);
semilogx(w, magnitude(W_closed_final, w), 'r');
grid on; ylabel('L, дБ');
subplot(2,1,2);
semilogx(w, phase(W_closed_final, w), 'r');
grid on; xlabel('\omega, рад/с'); ylabel('\phi, град');

function [W_ku, T1, T2, K, alpha, phi_max, L_base, L_lead] = lead(W, wc, gamma_req)
    Wc = evalfr(W, 1i*wc);
    gamma_base = 180 + rad2deg(angle(Wc));
    phi_max = gamma_req - gamma_base;
    alpha = (1 + sind(phi_max))/(1 - sind(phi_max));
    T2 = 1/(sqrt(alpha)*wc);
    T1 = alpha*T2;
    L_base = 20*log10(abs(Wc));
    L_lead = 10*log10(alpha);
    K = 10^(-(L_base + L_lead)/20);
    W_ku = K*tf([T1 1], [T2 1]);
end

function newFigure(name, height)
    fig = figure('Name', name, 'Color', 'w', 'Position', [100 100 850 height]);
    theme(fig, 'light');
    set(fig, 'defaultLineLineWidth', 1.5);
end

function L = magnitude(W, w)
    L = 20*log10(squeeze(bode(W, w)));
end

function phi = phase(W, w)
    [~, phi] = bode(W, w);
    phi = squeeze(phi);
end
