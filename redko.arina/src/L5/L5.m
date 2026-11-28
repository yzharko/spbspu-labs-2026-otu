clear;
clc;
close all;

%  W(s) = 32 / ((0.05s + 1)(0.1s + 1)(0.2s + 1))
num = 32;
den = conv(conv([0.05 1], [0.1 1]), [0.2 1]);
W_open = tf(num, den);
W_closed = feedback(W_open, 1);

char_eq = den + [0 0 0 num];
a3 = char_eq(1);
a2 = char_eq(2);
a1 = char_eq(3);
a0 = char_eq(4);

D1 = a2;
D2 = a2*a1 - a3*a0;
D3 = a0 * D2;

fprintf('Определители Гурвица: D1=%.4f  D2=%.4f  D3=%.4f\n', D1, D2, D3);

[Am, Pm, Wcg, Wcp] = margin(W_open);
fprintf('Исходная САУ: Pm = %.1f град на %.2f рад/с, Am = %.2f дБ на %.2f рад/с\n', ...
    Pm, Wcp, 20*log10(Am), Wcg);

figure('Name','Исходная САУ');
subplot(2,1,1);
step(W_closed, 1.2);
title('Переходная функция исходной САУ');
grid on;
subplot(2,1,2);
bode(W_open); grid on; title('ЛАЧХ/ЛФЧХ исходной разомкнутой САУ');

M_req = 1.19;
Tpp_req = 0.6;
sigma_req = 20;

gamma_req = asind(1/M_req);

Tp = [0.05 0.1 0.2];

% Нуль ПКУ компенсирует наибольшую постоянную времени объекта,
% полюс T2 выбирается в 5 раз меньше наименьшей постоянной объекта
% и выносит лишнее запаздывание за область частоты среза
T1 = max(Tp);
T2 = min(Tp) / 5;

% K выбирается из вариантов, удовлетворяющих M, sigma и Tпп,
% у которого установившееся значение h(inf) ближе к 1
w = logspace(-1, 3, 1000);
K = NaN;
h_best = 0;
for k = 0.001:0.001:0.2
    W_try = feedback(k * tf([T1 1], [T2 1]) * W_open, 1);
    if ~isstable(W_try)
        continue;
    end
    mag_try = squeeze(bode(W_try, w));
    M_try = max(mag_try) / dcgain(W_try);
    info_try = stepinfo(W_try);
    ok = M_try <= M_req && info_try.SettlingTime <= Tpp_req && info_try.Overshoot <= sigma_req;
    if ok && dcgain(W_try) > h_best
        h_best = dcgain(W_try);
        K = k;
    end
end
if isnan(K)
    error('При заданных Tp требования выполнить не удалось');
end

W_ku = K * tf([T1 1], [T2 1]);
fprintf('ПКУ: K = %.3f, T1 = %.3f с, T2 = %.3f с, h(inf) = %.3f\n', K, T1, T2, h_best);

asy = @(T) 20*log10(max(1, w*T));
L_isx = 20*log10(num) * ones(size(w));
for T = Tp
    L_isx = L_isx - asy(T);
end
L_zh = 20*log10(K*num) * ones(size(w));
for T = [Tp(Tp ~= T1) T2]
    L_zh = L_zh - asy(T);
end
L_ku = L_zh - L_isx;

W_final = W_ku * W_open;
W_closed_final = feedback(W_final, 1);

[Am2, Pm2, Wcg2, Wcp2] = margin(W_final);
fprintf('Скорректированная САУ: Pm = %.1f град (требуется >= %.1f), wc = %.2f рад/с\n', ...
    Pm2, gamma_req, Wcp2);

mag = squeeze(bode(W_closed_final, w));
M_fact = max(mag) / dcgain(W_closed_final);
fprintf('Фактическое M = %.4f (требуется <= %.2f)\n', M_fact, M_req);

info = stepinfo(W_closed_final);
Tpp_fact = info.SettlingTime;
sigma_fact = info.Overshoot;

yesNo = @(b) subsref({'НЕТ','ДА'}, struct('type','{}','subs',{{b+1}}));

fprintf('Соответствие: Tпп %.3f <= %.1f ? %s    sigma %.1f <= %.0f ? %s\n', ...
    Tpp_fact, Tpp_req, yesNo(Tpp_fact<=Tpp_req), sigma_fact, sigma_req, yesNo(sigma_fact<=sigma_req));

figure('Name','Асимптотические ЛАЧХ');
semilogx(w, L_isx, w, L_zh, w, L_ku, 'LineWidth', 1.5); grid on;
legend('Исходная L_{исх}','Желаемая L_{ж}','ПКУ L_{ку}','Location','best');
xlabel('\omega, рад/с'); ylabel('L, дБ');
title('Асимптотические ЛАЧХ — исходная, желаемая, ПКУ');

figure('Name','ЛАЧХ/ЛФЧХ систем');
bode(W_open, W_ku, W_final); grid on;
legend('Исходная W','ПКУ W_{ku}','Скорректированная W_{скор}','Location','best');
title('ЛАЧХ/ЛФЧХ — исходная, ПКУ, скорректированная');

figure('Name','Переходные характеристики');
step(W_closed, W_closed_final); grid on;
legend('Исходная САУ','Скорректированная САУ','Location','best');
title('Переходные характеристики замкнутых САУ');
