clear;
clc;
close all;

params.K = 3;
params.T = 1;
params.T1 = 1.2;
params.T2 = 1.5;
params.T1a = 3.5;

s = tf('s');
w = logspace(-2, 2, 1000);
w0 = 1 / sqrt(params.T);
w = w(abs(w - w0) > 1e-8);

links = [
    struct('n', 1, ...
           'name', 'Интегрирующее (идеальное)', ...
           'W', params.K/s)

    struct('n', 2, ...
           'name', 'Дифференцирующее (идеальное)', ...
           'W', params.K*s)

    struct('n', 3, ...
           'name', 'Усилительное (безынерционное)', ...
           'W', tf(params.K))

    struct('n', 4, ...
           'name', 'Апериодическое 1-го порядка', ...
           'W', params.K/(params.T*s + 1))

    struct('n', 5, ...
           'name', 'Апериодическое 2-го порядка', ...
           'W', params.K/(params.T2^2*s^2 + params.T1a*s + 1))

    struct('n', 6, ...
           'name', 'Колебательное', ...
           'W', params.K/(params.T2^2*s^2 + params.T1*s + 1))

    struct('n', 7, ...
           'name', 'Консервативное', ...
           'W', params.K/(params.T*s^2 + 1))

    struct('n', 8, ...
           'name', 'Инерционное интегрирующее', ...
           'W', params.K/(s*(params.T*s + 1)))

    struct('n', 9, ...
           'name', 'Реальное дифференцирующее', ...
           'W', params.K*s/(params.T*s + 1))

    struct('n', 10, ...
           'name', 'Форсирующее 1-го порядка', ...
           'W', params.K*(params.T*s + 1))

    struct('n', 11, ...
           'name', 'Изодромное', ...
           'W', params.K*(params.T*s + 1)/s)
];

%% Построение АФЧХ динамических звеньев

figure('Name', 'АФЧХ динамических звеньев', 'Position', [0 0 1200 900]);
tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(links)
    nexttile;

    W = links(i).W;
    n = links(i).n;

    if n == 3
        plot(params.K, 0, 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b');
        grid on;
    elseif n == 7
        Re = params.K ./ (1 - params.T * w.^2);
        plot(Re, zeros(size(Re)), 'b', 'LineWidth', 1);
        grid on;
    else
        nyquist(W, w);
    end

    title(sprintf('№%d. %s', n, links(i).name), 'FontSize', 12);

    xlabel('Re(W)');
    ylabel('Im(W)');
end
sgtitle('АФЧХ динамических звеньев', 'FontSize', 14, 'FontWeight', 'bold');

for i = 1:length(links)
    [mag, phase] = bode(links(i).W, w);
    
    mag = squeeze(mag);
    phase = squeeze(phase);
    
    links(i).L = 20 * log10(mag);
    links(i).phase = phase;
end

%% Построение ЛАЧХ динамических звеньев

figure('Name', 'ЛАЧХ динамических звеньев', 'Position', [0 0 1200 900]);
tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(links)
    nexttile;
    
    semilogx(w, links(i).L, 'LineWidth', 1);
    grid on;

    title(sprintf('№%d. %s', links(i).n, links(i).name), 'FontSize', 12);
    xlabel('Частота, рад/с');
    ylabel('Усиление звена, дБ');
end
sgtitle('ЛАЧХ динамических звеньев', 'FontSize', 14, 'FontWeight', 'bold');

%% Построение ЛФЧХ динамических звеньев

figure('Name', 'ЛФЧХ динамических звеньев', 'Position', [0 0 1200 900]);
tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(links)
    nexttile;
    
    semilogx(w, links(i).phase, 'LineWidth', 1);
    grid on;

    title(sprintf('№%d. %s', links(i).n, links(i).name), 'FontSize', 12);
    xlabel('Частота, рад/с');
    ylabel('Фазовый сдвиг, градусы');
end
sgtitle('ЛФЧХ динамических звеньев', 'FontSize', 14, 'FontWeight', 'bold');

%% Исследование влияния K на ЛАЧХ

K_mult = [0.5 1 1.5 2 5];
labelsK = cell(1, length(K_mult));

for i = 1:length(K_mult)
    labelsK{i} = sprintf('K*%.1f (K=%.1f)', K_mult(i), params.K*K_mult(i));
end

figure('Name', 'Влияние коэффициента усиления на ЛАЧХ', 'Position', [0 0 1200 900]);
for n = 1:length(links)
    subplot(3, 4, n);
    hold on;

    for i = 1:length(K_mult)
        W_new = links(n).W * K_mult(i); 
        
        [mag, ~] = bode(W_new, w);
        semilogx(w, 20*log10(squeeze(mag)), 'LineWidth', 1);
    end

    hold off;
    grid on;
    set(gca, 'XScale', 'log');
    xlabel('Частота, рад/с');
    ylabel('Уровень, дБ');
    title(sprintf('№%d. %s', n, links(n).name), 'FontSize', 9);
end
subplot(3, 4, 12);
add_legend_plot(labelsK);
sgtitle('Влияние коэффициента усиления K на ЛАЧХ');

%% Исследование влияния K на ЛФЧХ

figure('Name', 'Влияние коэффициента усиления на ЛФЧХ', 'Position', [0 0 1200 900]);
for n = 1:length(links)
    subplot(3, 4, n);
    hold on;

    for i = 1:length(K_mult)
        W_new = links(n).W * K_mult(i);
        [~, phase] = bode(W_new, w);
        semilogx(w, squeeze(phase), 'LineWidth', 1);
    end

    hold off;
    grid on;
    set(gca, 'XScale', 'log');
    xlabel('Частота, рад/с');
    ylabel('Фазовый сдвиг, град');
    title(sprintf('№%d. %s', n, links(n).name), 'FontSize', 9);
end
subplot(3, 4, 12);
add_legend_plot(labelsK);
sgtitle('Влияние коэффициента усиления K на ЛФЧХ');

%% Исследование влияния T на ЛАЧХ

T_mult = [0.1 0.5 1 2 5]; 
T0 = params.T;
labelsT = cell(1, length(T_mult));
for i = 1:length(T_mult)
    labelsT{i} = sprintf('T*%.1f (T=%.1f)', T_mult(i), T_mult(i)*params.T);
end

WT = cell(length(links), length(T_mult));

for i = 1:length(T_mult)
    a = T_mult(i);
    T_new = a * T0;
    T1_new = a * params.T1;
    T2_new = a * params.T2;
    T1a_new = a * params.T1a;
    
    WT{1, i}  = params.K / s;
    WT{2, i}  = params.K * s;
    WT{3, i}  = tf(params.K);
    WT{4, i}  = params.K / (T_new*s + 1);
    WT{5, i}  = params.K / (T2_new^2*s^2 + T1a_new*s + 1);
    WT{6, i}  = params.K / (T2_new^2*s^2 + T1_new*s + 1);
    WT{7, i}  = params.K / (T_new*s^2 + 1);
    WT{8, i}  = params.K / (s*(T_new*s + 1));
    WT{9, i}  = params.K*s / (T_new*s + 1);
    WT{10, i} = params.K * (T_new*s + 1);
    WT{11, i} = params.K * (T_new*s + 1) / s;
end

figure('Name', 'Влияние постоянных времени на ЛАЧХ', 'Position', [0 0 1200 900]);
for n = 1:length(links)
    subplot(3, 4, n);
    hold on;

    for i = 1:length(T_mult)
        [mag, ~] = bode(WT{n, i}, w);
        semilogx(w, 20*log10(squeeze(mag)), 'LineWidth', 1);
    end

    hold off;
    grid on;
    set(gca, 'XScale', 'log');
    xlabel('Частота, рад/с');
    ylabel('Уровень, дБ');
    title(sprintf('№%d. %s', n, links(n).name), 'FontSize', 9);
end
subplot(3, 4, 12);
add_legend_plot(labelsT);
sgtitle('Влияние постоянных времени T на ЛАЧХ');

%% Исследование влияния T на ЛФЧХ

figure('Name', 'Влияние постоянных времени на ЛФЧХ', 'Position', [0 0 1200 900]);
for n = 1:length(links)
    subplot(3, 4, n);
    hold on;

    for i = 1:length(T_mult)
        [~, phase] = bode(WT{n, i}, w);
        semilogx(w, squeeze(phase), 'LineWidth', 1);
    end

    hold off;
    grid on;
    set(gca, 'XScale', 'log');
    xlabel('Частота, рад/с');
    ylabel('Фазовый сдвиг, град');
    title(sprintf('№%d. %s', n, links(n).name), 'FontSize', 9);
end
subplot(3, 4, 12);
add_legend_plot(labelsT);
sgtitle('Влияние постоянных времени T на ЛФЧХ');

function add_legend_plot(labels)
    hold on;
    for i = 1:length(labels)
        plot(nan, nan, 'LineWidth', 1);
    end
    axis off;
    legend(labels, 'Location', 'northwest', 'FontSize', 10);
end
