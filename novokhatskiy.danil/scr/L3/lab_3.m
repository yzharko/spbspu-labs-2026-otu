clc;
clear;
close all;

s = tf('s');

K0 = 1.1;
T0 = 0.2;
T1_0 = 0.3;
T2_0 = 0.5;
T1_ap = 1.5;

omega = logspace(-2, 2, 2000);
omegaN = logspace(log10(0.05), log10(20), 4000);

names = {
    'Интегрирующее';
    'Дифференцирующее';
    'Усилительное';
    'Апериодическое 1-го';
    'Апериодическое 2-го';
    'Колебательное';
    'Консервативное';
    'Инерционное интегр.';
    'Реальное дифф.';
    'Форсирующее';
    'Изодромное'};

W = cell(1, 11);
W{1}  = K0/s;
W{2}  = K0*s;
W{3}  = tf(K0);
W{4}  = K0/(T0*s + 1);
W{5}  = K0/(T2_0^2*s^2 + T1_ap*s + 1);
W{6}  = K0/(T2_0^2*s^2 + T1_0*s + 1);
W{7}  = K0/(T0*s^2 + 1); 
W{8}  = K0/(s*(T0*s + 1));
W{9}  = K0*s/(T0*s + 1);
W{10} = K0*(T0*s + 1);
W{11} = K0*(T0*s + 1)/s;

figure('Name', 'АФЧХ всех звеньев', 'Position', [50 50 1450 850]);
for n = 1:length(W)
    [re, im] = nyquist(W{n}, omegaN);
    re = squeeze(re);
    im = squeeze(im);
    if n == 1 || n == 2
        re(:) = 0; 
    elseif n == 7
        re = K0./(1 - T0*omegaN(:).^2);
        im = zeros(size(re));
        j = find(omegaN > 1/sqrt(T0), 1, 'first');
        re(j) = NaN; 
    elseif n == 10
        re(:) = K0;
    elseif n == 11
        re(:) = K0*T0;  
        im = -K0./omegaN(:);
    end
    subplot(3, 4, n);
    if n == 3
        plot(K0, 0, 'bo', 'MarkerFaceColor', 'b');
    else
        plot(re, im, 'b', re, -im, 'b');
    end
    grid on;
    xlabel('Вещественная часть');
    ylabel('Мнимая часть');
    title(sprintf('%d. %s', n, names{n}), 'FontSize', 9);
end
sgtitle('АФЧХ типовых динамических звеньев');

figure('Name', 'ЛАЧХ всех звеньев', 'Position', [70 70 1450 850]);
for n = 1:length(W)
    [mag, ~] = bode(W{n}, omega);
    subplot(3, 4, n);
    semilogx(omega, 20*log10(squeeze(mag)), 'b');
    grid on;
    xlabel('Частота, рад/с');
    ylabel('Уровень, дБ');
    title(sprintf('%d. %s', n, names{n}), 'FontSize', 9);
end
sgtitle('ЛАЧХ типовых динамических звеньев');

figure('Name', 'ЛФЧХ всех звеньев', 'Position', [90 90 1450 850]);
for n = 1:length(W)
    [~, phase] = bode(W{n}, omega);
    subplot(3, 4, n);
    semilogx(omega, squeeze(phase), 'm');
    grid on;
    xlabel('Частота, рад/с');
    ylabel('Фазовый сдвиг, град');
    title(sprintf('%d. %s', n, names{n}), 'FontSize', 9);
end
sgtitle('ЛФЧХ типовых динамических звеньев');

K_values = [0.5 K0 2 3 5];
labelsK = cell(1, length(K_values));
for i = 1:length(K_values)
    labelsK{i} = sprintf('Коэффициент усиления = %.2f', K_values(i));
end

figure('Name', 'Влияние коэффициента усиления на ЛАЧХ', 'Position', [110 110 1450 850]);
for n = 1:length(W)
    WK = cell(1, length(K_values));
    subplot(3, 4, n);
    hold on;
    for i = 1:length(K_values)
        WK{i} = (K_values(i)/K0)*W{n};
        [mag, ~] = bode(WK{i}, omega);
        semilogx(omega, 20*log10(squeeze(mag)));
    end
    hold off;
    set(gca, 'XScale', 'log');
    grid on;
    xlabel('Частота, рад/с');
    ylabel('Уровень, дБ');
    title(sprintf('%d. %s', n, names{n}), 'FontSize', 9);
end
subplot(3, 4, 12);
hold on;
for i = 1:length(K_values)
    plot(nan, nan);
end
axis off;
legend(labelsK, 'Location', 'northwest');
sgtitle('Влияние коэффициента усиления на ЛАЧХ');

figure('Name', 'Влияние коэффициента усиления на ЛФЧХ', 'Position', [130 130 1450 850]);
for n = 1:length(W)
    subplot(3, 4, n);
    hold on;
    for i = 1:length(K_values)
        WK{i} = (K_values(i)/K0)*W{n};
        [~, phase] = bode(WK{i}, omega);
        semilogx(omega, squeeze(phase));
    end
    hold off;
    set(gca, 'XScale', 'log');
    grid on;
    xlabel('Частота, рад/с');
    ylabel('Фазовый сдвиг, град');
    title(sprintf('%d. %s', n, names{n}), 'FontSize', 9);
end
subplot(3, 4, 12);
hold on;
for i = 1:length(K_values)
    plot(nan, nan);
end
axis off;
legend(labelsK, 'Location', 'northwest');
sgtitle('Влияние коэффициента усиления на ЛФЧХ');

T_factors = [0.5 1 2 3 5];
T_values = T0*T_factors;
WT = cell(length(W), length(T_factors));
labelsT = cell(1, length(T_factors));

for i = 1:length(T_factors)
    a = T_factors(i);
    T = a*T0;
    T1 = a*T1_0;
    T2 = a*T2_0;
    T1a = a*T1_ap;
    WT{1, i}  = K0/s;
    WT{2, i}  = K0*s;
    WT{3, i}  = tf(K0);
    WT{4, i}  = K0/(T*s + 1);
    WT{5, i}  = K0/(T2^2*s^2 + T1a*s + 1);
    WT{6, i}  = K0/(T2^2*s^2 + T1*s + 1);
    WT{7, i}  = K0/(T*s^2 + 1);
    WT{8, i}  = K0/(s*(T*s + 1));
    WT{9, i}  = K0*s/(T*s + 1);
    WT{10, i} = K0*(T*s + 1);
    WT{11, i} = K0*(T*s + 1)/s;
    labelsT{i} = sprintf('Масштаб времени = %.1f', a);
end

figure('Name', 'Влияние постоянных времени на ЛАЧХ', 'Position', [150 150 1450 850]);
for n = 1:length(W)
    subplot(3, 4, n);
    hold on;
    for i = 1:length(T_factors)
        [mag, ~] = bode(WT{n, i}, omega);
        semilogx(omega, 20*log10(squeeze(mag)));
    end
    hold off;
    set(gca, 'XScale', 'log');
    grid on;
    xlabel('Частота, рад/с');
    ylabel('Уровень, дБ');
    title(sprintf('%d. %s', n, names{n}), 'FontSize', 9);
end
subplot(3, 4, 12);
hold on;
for i = 1:length(T_factors)
    plot(nan, nan);
end
axis off;
legend(labelsT, 'Location', 'northwest');
sgtitle('Влияние постоянных времени на ЛАЧХ');

figure('Name', 'Влияние постоянных времени на ЛФЧХ', 'Position', [170 170 1450 850]);
for n = 1:length(W)
    subplot(3, 4, n);
    hold on;
    for i = 1:length(T_factors)
        [~, phase] = bode(WT{n, i}, omega);
        semilogx(omega, squeeze(phase));
    end
    hold off;
    set(gca, 'XScale', 'log');
    grid on;
    xlabel('Частота, рад/с');
    ylabel('Фазовый сдвиг, град');
    title(sprintf('%d. %s', n, names{n}), 'FontSize', 9);
end
subplot(3, 4, 12);
hold on;
for i = 1:length(T_factors)
    plot(nan, nan);
end
axis off;
legend(labelsT, 'Location', 'northwest');
sgtitle('Влияние постоянных времени на ЛФЧХ');

figure('Name', 'Влияние параметров на АФЧХ апериодического звена', ...
    'Position', [180 180 1200 500]);
subplot(1, 2, 1);
hold on;
for i = 1:length(K_values)
    WK_ap = K_values(i)/(T0*s + 1);
    [re, im] = nyquist(WK_ap, omegaN);
    h = plot(squeeze(re), squeeze(im));
    plot(squeeze(re), -squeeze(im), 'Color', h.Color, ...
        'HandleVisibility', 'off');
end
hold off;
grid on;
xlabel('Вещественная часть'); ylabel('Мнимая часть');
title('Изменение коэффициента усиления');
legend(labelsK, 'Location', 'best');

subplot(1, 2, 2);
hold on;
for i = 1:length(T_values)
    WT_ap = K0/(T_values(i)*s + 1);
    [re, im] = nyquist(WT_ap, omegaN);
    h = plot(squeeze(re), squeeze(im));
    plot(squeeze(re), -squeeze(im), 'Color', h.Color, ...
        'HandleVisibility', 'off');
end
hold off;
grid on;
xlabel('Вещественная часть'); ylabel('Мнимая часть');
title('Изменение постоянной времени');
legend(labelsT, 'Location', 'best');
sgtitle('Сравнение АФЧХ апериодического звена 1-го порядка');

T1_ap_values = [1.0 T1_ap 2.0];
T2_ap_values = [0.3 0.4 T2_0];
T1_osc_values = [0.2 T1_0 0.5];
T2_osc_values = [0.4 T2_0 0.6];

for n = [5 6]
    if n == 5
        values1 = T1_ap_values;
        values2 = T2_ap_values;
        T1_base = T1_ap;
    else
        values1 = T1_osc_values;
        values2 = T2_osc_values;
        T1_base = T1_0;
    end

    W1 = cell(1, length(values1));
    W2 = cell(1, length(values2));
    labels1 = cell(1, length(values1));
    labels2 = cell(1, length(values2));
    for i = 1:length(values1)
        W1{i} = K0/(T2_0^2*s^2 + values1(i)*s + 1);
        W2{i} = K0/(values2(i)^2*s^2 + T1_base*s + 1);
        labels1{i} = sprintf('Первая постоянная = %.2f с', values1(i));
        labels2{i} = sprintf('Вторая постоянная = %.2f с', values2(i));
    end

    figure('Name', ['Влияние постоянных времени: ' names{n}], ...
        'Position', [190 190 1200 800]);
    for i = 1:length(values1)
        [mag, phase] = bode(W1{i}, omega);
        subplot(2, 2, 1);
        hold on;
        semilogx(omega, 20*log10(squeeze(mag)));
        subplot(2, 2, 2);
        hold on;
        semilogx(omega, squeeze(phase));
    end
    for i = 1:length(values2)
        [mag, phase] = bode(W2{i}, omega);
        subplot(2, 2, 3);
        hold on;
        semilogx(omega, 20*log10(squeeze(mag)));
        subplot(2, 2, 4);
        hold on;
        semilogx(omega, squeeze(phase));
    end
    for k = 1:4
        subplot(2, 2, k);
        hold off;
        set(gca, 'XScale', 'log');
        grid on;
        xlabel('Частота, рад/с');
        if k == 1 || k == 3
            ylabel('Уровень, дБ');
        else
            ylabel('Фазовый сдвиг, град');
        end
        if k <= 2
            legend(labels1, 'Location', 'best');
        else
            legend(labels2, 'Location', 'best');
        end
    end
    subplot(2, 2, 1); title('Первая постоянная времени: ЛАЧХ');
    subplot(2, 2, 2); title('Первая постоянная времени: ЛФЧХ');
    subplot(2, 2, 3); title('Вторая постоянная времени: ЛАЧХ');
    subplot(2, 2, 4); title('Вторая постоянная времени: ЛФЧХ');
    sgtitle(['Звено № ' num2str(n) ': ' names{n}]);
end
