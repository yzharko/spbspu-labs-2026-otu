clc;
clear;
close all;

s = tf('s');

K0 = 0.9;
T0 = 0.6;
T1 = 0.8;
T2 = 0.7;
T1_aperiodic = 1.5;

W_names = {"Интегрирующее", ...
    "Дифференцирующее", ...
    "Усилительное", ...
    "Апериодическое 1-го порядка", ...
    "Апериодическое 2-го порядка", ...
    "Колебательное", ...
    "Консервативное", ...
    "Инерционное интегрирующее", ...
    "Реальное дифференцирующее", ...
    "Форсирующее 1-го порядка", ...
    "Изодромное"};

for i = 1:numel(W_names)
    % Построение ЛАЧХ, ЛФЧХ и АФЧХ
    W = W_by_index(i, K0, T0, T1, T2, T1_aperiodic, s);
    if isnumeric(W)
        W = tf(W);
    end
    figure;
    bode(W);
    grid on;
    title(['ЛАЧХ и ЛФЧХ: ' W_names{i}]);
end

for i = 1:numel(W_names)
    W = W_by_index(i, K0, T0, T1, T2, T1_aperiodic, s);
    if isnumeric(W)
        W = tf(W);
    end
    figure;
    nyquist(W);
    grid on;
    title(['АФЧХ: ' W_names{i}]);
end

% Определяем влияния коэффициентов
K_values = [0.5 1 2 3 5];
T_values = [0.1 0.5 1 2 5];

% T фиксировано, меняем K
for j = 1:numel(W_names)
    WK = cell(1, length(K_values));
    labelsK = cell(1, length(K_values));

    for i = 1:length(K_values)
        K = K_values(i);
        WK{i} = W_by_index(j, K, T0, T1, T2, T1_aperiodic, s);
        if isnumeric(WK{i})
            WK{i} = tf(WK{i});
        end
        labelsK{i} = sprintf('K = %.2g', K);
    end
    figure;
    nyquist(WK{:});
    grid on;
    legend(labelsK, 'Location', 'best');
    title(['АФЧХ: ' W_names{j} ' — влияние K (T фиксировано)']);

    figure;
    bode(WK{:});
    grid on;
    legend(labelsK, 'Location', 'best');
    title(['ЛАЧХ и ЛФЧХ: ' W_names{j} ' — влияние K (T фиксировано)']);
end

% K фиксировано, меняем T
for j = 1:numel(W_names)
    WT = cell(1, length(T_values));
    labelsT = cell(1, length(T_values));

    for i = 1:length(T_values)
        T = T_values(i);
        WT{i} = W_by_index(j, K0, T, T1, T2, T1_aperiodic, s);
        if isnumeric(WT{i})
            WT{i} = tf(WT{i});
        end
        labelsT{i} = sprintf('T = %.2g', T);
    end
    figure;
    nyquist(WT{:});
    grid on;
    legend(labelsT, 'Location', 'best');
    title(['АФЧХ: ' W_names{j} ' — влияние T (K фиксировано)']);

    figure;
    bode(WT{:});
    grid on;
    legend(labelsT, 'Location', 'best');
    title(['ЛАЧХ и ЛФЧХ: ' W_names{j} ' — влияние T (K фиксировано)']);
end

% Исследование влияния T1 и T2 для звеньев второго порядка
index = [5 6];

T1_values = {
    [1.4 1.5 2 3 5];      % Для апериодического
    [0.1 0.3 0.5 0.8 1.3] % Для колебательного
    };

T2_values = {
    [0.1 0.3 0.5 0.7];    % Для апериодического
    [0.41 0.5 0.7 1 1.5]  % Для колебательного
    };

for j = 1:length(index)

    idx = index(j);
    name = W_names{idx};

    values = T1_values{j};

    W_T1 = cell(1, length(values));
    labels_T1 = cell(1, length(values));

    for n = 1:length(values)
        T1x = values(n);
        W_T1{n} = W_by_index(idx, K0, T0, T1x, T2, T1x, s);
        W_T1{n} = tf(W_T1{n});

        labels_T1{n} = sprintf('T1 = %.2g', T1x);
    end

    figure;
    bode(W_T1{:});
    grid on;
    legend(labels_T1, 'Location', 'best');
    title(['ЛАЧХ и ЛФЧХ: ' name ' — влияние T1']);

    figure;
    nyquist(W_T1{:});
    grid on;
    legend(labels_T1, 'Location', 'best');
    title(['АФЧХ: ' name ' — влияние T1']);

    values = T2_values{j};

    W_T2 = cell(1, length(values));
    labels_T2 = cell(1, length(values));

    for n = 1:length(values)
        T2x = values(n);
        W_T2{n} = W_by_index(idx, K0, T0, T1, T2x, T1_aperiodic, s);
        W_T2{n} = tf(W_T2{n});
        labels_T2{n} = sprintf('T2 = %.2g', T2x);
    end

    figure;
    bode(W_T2{:});
    grid on;
    legend(labels_T2, 'Location', 'best');
    title(['ЛАЧХ и ЛФЧХ: ' name ' — влияние T2']);

    figure;
    nyquist(W_T2{:});
    grid on;
    legend(labels_T2, 'Location', 'best');
    title(['АФЧХ: ' name ' — влияние T2']);

end

% Небольшая вспомогательная аналитика
fprintf('\nСопрягающие частоты при K = %.2f:\n', K0);

for i = 1:numel(W_names)
    fprintf('\nЗвено %d:\n', i);
    if i == 1 || i == 2 || i == 3
        continue;
    elseif i == 4 || i == 8 || i == 9 || i == 10 || i == 11
        for n = 1:length(T_values)
            T = T_values(n);
            omega_c = 1/T;
            fprintf('T = %.2f c, omega_c = %.3f рад/с\n', T, omega_c);
        end

    elseif i == 5 || i == 6
        omega_n = 1/T2;
        fprintf('T = %.2f c, omega_n = %.3f рад/с\n', T, omega_n);

    elseif i == 7
        for n = 1:length(T_values)
            T = T_values(n);
            omega_0 = 1/sqrt(T);
            fprintf('T = %.2f c, omega_0 = %.3f рад/с\n', T, omega_0);
        end
    end
end

function Wi = W_by_index(i, Kx, Tx, T1x, T2x, T1ax, sx)
if i == 1
    Wi = Kx / sx;
elseif i == 2
    Wi = Kx * sx;
elseif i == 3
    Wi = Kx;
elseif i == 4
    Wi = Kx / (Tx*sx + 1);
elseif i == 5
    Wi = Kx / (T2x^2*sx^2 + T1ax*sx + 1);
elseif i == 6
    Wi = Kx / (T2x^2*sx^2 + T1x*sx + 1);
elseif i == 7
    Wi = Kx / (Tx*sx^2 + 1);
elseif i == 8
    Wi = Kx / (sx*(Tx*sx + 1));
elseif i == 9
    Wi = Kx*sx / (Tx*sx + 1);
elseif i == 10
    Wi = Kx * (Tx*sx + 1);
elseif i == 11
    Wi = Kx * (Tx*sx + 1) / sx;
end
end

