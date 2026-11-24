clear; 
clc; 
close all;

%% 1. Исследование исходной САУ
%  W(s) = 228 / 0.00014*s^4 + 0.0088*s^3 + 0.28*s^2 + s
%  Tпп = 0.7 с, M = 1.3, sigma = 10 %
num = 228;
den = [0.00014  0.0088  0.28  2  0];
W_open = tf(num, den);
W_closed = feedback(W_open, 1);

% Проверка устойчивости по Гурвицу
char_eq = den + [0 0 0 0 num];
a4 = char_eq(5); 
a3 = char_eq(4); 
a2 = char_eq(3); 
a1 = char_eq(2); 
a0 = char_eq(1);

D1 = a3;
D2 = a3*a2 - a4*a1;
D3 = a3*a2*a1 - a3^2*a0 - a4*a1^2;
D4 = a0 * D3;

fprintf('Определители Гурвица: D1=%.2f  D2=%.2f  D3=%.2f  D4=%.2f\n', D1, D2, D3, D4);

% Исходные запасы устойчивости
[Am, Pm, Wcg, Wcp] = margin(W_open);

figure('Name','Исходная САУ');
subplot(2,1,1);
step(W_closed); 
title('Переходная функция исходной САУ');
grid on;
subplot(2,1,2);
bode(W_open); grid on; title('ЛАЧХ/ЛФЧХ исходной разомкнутой САУ');

%% 2. Расчёт параметров желаемой ЛАЧХ
M_req = 1.3;
Tpp_req = 0.7;
sigma_req = 20;

% Требуемый запас по фазе из M
gamma_req = asind(1/M_req);

% Эмпирическим путём воздействуем на wc и убеждаемся, 
% что результат попадает в заданные требования
wc = 10;

% Исходная фаза на частоте wc
[mag_i, ph_i] = bode(W_open, wc);
mag_i = squeeze(mag_i); 
ph_i = squeeze(ph_i);
gamma_base = 180 + ph_i;

%% 3. Рассчёт ПКУ
% W_ku(s) = K * (T1*s + 1) / (T2*s + 1)

% Требуемое добавочное фазовое опережение
delta_phi = gamma_req - gamma_base;
alpha = (1 + sind(delta_phi)) / (1 - sind(delta_phi));

% w_m = 1/sqrt(T1*T2) = wc   =>   T1*T2 = 1/wc^2,  T1 = alpha*T2
T2 = 1 / (wc * sqrt(alpha));
T1 = alpha * T2;

% Определяем K из условия L_желаемое(w) = 0
L_wc = 20*log10(mag_i);
% Амплитудный вклад звена (T1 s+1)/(T2 s+1) на wc
L_zv_wc = 20*log10(sqrt(1+(T1*wc)^2) / sqrt(1+(T2*wc)^2));
L_K = 0 - L_wc - L_zv_wc;
K = 10^(L_K/20);

% Передаточная функция КУ
W_ku = K * tf([T1 1], [T2 1]);

%% 4. Скорректированная САУ
W_final = W_ku * W_open;
W_closed_final = feedback(W_final, 1);

% Запасы устойчивости скорректированной системы
[Am2, Pm2, Wcg2, Wcp2] = margin(W_final);

M_fact = 1 / sind(Pm2);
fprintf('Фактическое M = 1/sin(Pm) = %.2f (требуется <= %.2f)\n', M_fact, M_req);

%% 5. Показатели качества переходного процесса
info = stepinfo(W_closed_final);
Tpp_fact = info.SettlingTime;
sigma_fact = info.Overshoot;

% Локальная функция для форматирования булевых значений
yesNo = @(b) subsref({'НЕТ','ДА'}, struct('type','{}','subs',{{b+1}}));

fprintf('Соответствие: Tпп %.3f <= %.1f ? %s    sigma %.1f <= %.0f ? %s\n', ...
    Tpp_fact, Tpp_req, yesNo(Tpp_fact<=Tpp_req), sigma_fact, sigma_req, yesNo(sigma_fact<=sigma_req));

%% 6. ЛАЧХ/ЛФЧХ исходной и скорректированной САУ
figure('Name','ЛАЧХ/ЛФЧХ систем');
subplot(2,1,1);
bode(W_open, W_ku, W_final); grid on;
legend('Исходная W','ПКУ W_{ku}','Скорректированная W_{скор}','Location','best');
title('ЛАЧХ (Боде) — исходная, ПКУ, скорректированная');
subplot(2,1,2);
bode(W_open, W_ku, W_final); grid on;
title('ЛФЧХ');

% Переходные характеристики до и после
figure('Name','Переходные характеристики');
step(W_closed, W_closed_final); grid on;
legend('Исходная САУ','Скорректированная САУ','Location','best');
title('Переходные характеристики замкнутых САУ');
