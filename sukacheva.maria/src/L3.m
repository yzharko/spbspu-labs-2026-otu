clc;
clear;
close all;

K_val = 1.6;
T = 0.8;
T1 = 1.0;
T1_2_aperiodic = 3.0;
T2 = 1.1;

s = tf('s');

names = {
    '1. Интегрирующее'
    '2. Дифференцирующее (идеальное)'
    '3. Усилительное (безынерционное)'
    '4. Апериодическое 1-го порядка'
    '5. Апериодическое 2-го порядка'
    '6. Колебательное'
    '7. Консервативное'
    '8. Инерционное интегрирующее'
    '9. Реальное дифференцирующее'
    '10. Форсирующее 1-го порядка'
    '11. Изодромное'
};

w = logspace(-2, 2, 1000);

Wc_all = {
    K_val / s
    K_val * s
    tf(K_val)
    K_val / (T * s + 1)
    K_val / (T2^2 * s^2 + T1_2_aperiodic * s + 1)
    K_val / (T2^2 * s^2 + T1 * s + 1)
    K_val / (T * s^2 + 1)
    K_val / (s * (T * s + 1))
    (K_val * s) / (T * s + 1)
    K_val * (T * s + 1)
    (K_val * (T * s + 1)) / s
};

figure('Position', [100 100 1300 900]);
for i = 1:numel(Wc_all)
    subplot(4, 3, i);
    [mag, ~, ~] = bode(Wc_all{i}, w);
    semilogx(w, 20*log10(squeeze(mag)), 'LineWidth', 1.5);
    grid on;
    title(names{i}, 'FontSize', 9);
    xlabel('\omega, рад/с');
    ylabel('L, дБ');
end
sgtitle('ЛАЧХ');

figure('Position', [100 100 1300 900]);
for i = 1:numel(Wc_all)
    subplot(4, 3, i);
    [~, phase, ~] = bode(Wc_all{i}, w);
    semilogx(w, squeeze(phase), 'LineWidth', 1.5);
    grid on;
    title(names{i}, 'FontSize', 9);
    xlabel('\omega, рад/с');
    ylabel('\psi, град');
end
sgtitle('ЛФЧХ');

figure('Units', 'normalized', 'Position', [0 0 1 1]);
for i = 1:numel(Wc_all)
    subplot(3, 4, i);
    nyquist(Wc_all{i});
    grid on;
    title(names{i}, 'FontSize', 9);
end
sgtitle('АФЧХ');

mult = [1/5 1/2 1 2 5];

Kv = [K_val * mult, K_val * ones(1, 5)];
mv = [ones(1, 5), mult];

Wv = cell(1, 10);
for j = 1:10
    K = Kv(j);
    m = mv(j);
    Wv{j} = {
        K / s
        K * s
        tf(K)
        K / (T*m * s + 1)
        K / ((T2*m)^2 * s^2 + T1_2_aperiodic*m * s + 1)
        K / ((T2*m)^2 * s^2 + T1*m * s + 1)
        K / (T*m * s^2 + 1)
        K / (s * (T*m * s + 1))
        (K * s) / (T*m * s + 1)
        K * (T*m * s + 1)
        (K * (T*m * s + 1)) / s
    };
end

labelsK = cell(1, 5);
labelsT = cell(1, 5);
for j = 1:5
    labelsK{j} = sprintf('K = %.2g', K_val * mult(j));
    labelsT{j} = sprintf('T = %.2g', T * mult(j));
end

plot_tiles(Wv(1:5),  names, w, 'mag',   labelsK, ['ЛАЧХ: T = ', num2str(T), ' (фикс.), меняется K']);
plot_tiles(Wv(1:5),  names, w, 'phase', labelsK, ['ЛФЧХ: T = ', num2str(T), ' (фикс.), меняется K']);
plot_tiles(Wv(6:10), names, w, 'mag',   labelsT, ['ЛАЧХ: K = ', num2str(K_val), ' (фикс.), меняется T']);
plot_tiles(Wv(6:10), names, w, 'phase', labelsT, ['ЛФЧХ: K = ', num2str(K_val), ' (фикс.), меняется T']);

function plot_tiles(Wset, names, w, kind, labels, ttl)

figure('Units', 'normalized', 'Position', [0 0 1 1]);

for i = 1:numel(names)

    subplot(4, 3, i);
    hold on;
    grid on;

    for j = 1:numel(Wset)
        [mag, phase, ~] = bode(Wset{j}{i}, w);
        if strcmp(kind, 'mag')
            y = 20*log10(squeeze(mag));
        else
            y = squeeze(phase);
        end
        semilogx(w, y, 'LineWidth', 1.3);
    end

    set(gca, 'XScale', 'log');
    xlim([w(1) w(end)]);

    title(names{i}, 'FontSize', 9);
    xlabel('\omega, рад/с');
    if strcmp(kind, 'mag')
        ylabel('L, дБ');
    else
        ylabel('\psi, град');
    end

    hold off;
end

subplot(4, 3, 12);
hold on;
for j = 1:numel(labels)
    plot(nan, nan, 'LineWidth', 1.5);
end
axis off;
legend(labels, 'Location', 'best');
hold off;

sgtitle(ttl);

end
