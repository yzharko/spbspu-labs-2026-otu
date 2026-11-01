clc;
clear;
close all;

s = tf('s');

K0 = 1;
T0 = 1;
W = K0/(T0*s + 1);

%% Построение ЛАЧХ, ЛФЧХ и АФЧХ
figure;
bode(W);
grid on;
title('ЛАЧХ и ЛФЧХ апериодического звена');

figure;
nyquist(W);
grid on;
title('АФЧХ апериодического звена');

%% Определяем влияния коэффициентов
K_values = [0.5 1 2 3 5];
T_values = [0.1 0.5 1 2 5];

% T фиксировано, меняем K
WK = cell(1, length(K_values));
labelsK = cell(1, length(K_values));

for i = 1:length(K_values)
    K = K_values(i);
    WK{i} = K/(T0*s + 1);
    labelsK{i} = sprintf('K = %.2g', K);
end

figure;
nyquist(WK{:});
grid on;
legend(labelsK, 'Location', 'best');

figure;
bode(WK{:});
grid on;
legend(labelsK, 'Location', 'best');

% K фиксировано, меняем T
WT = cell(1, length(T_values));
labelsT = cell(1, length(T_values));

for i = 1:length(T_values)
    T = T_values(i);
    WT{i} = K0/(T*s + 1);
    labelsT{i} = sprintf('T = %.2g', T);
end

figure;
nyquist(WT{:});
grid on;
legend(labelsT, 'Location', 'best');

figure;
bode(WT{:});
grid on;
legend(labelsT, 'Location', 'best');

%% Небольшая вспомогательная аналитика
fprintf('\nСопрягающие частоты при K = %.2f:\n', K0);
for i = 1:length(T_values)
    T = T_values(i);
    omega_c = 1/T;

    fprintf('T = %.2f c, omega_c = %.3f рад/с\n', ...
        T, omega_c);
end
