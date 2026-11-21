clear; clc; close all;

K  = 0.7;
T  = 0.4;
T1 = 0.5;
T2 = 0.7;
T1_aper = 1.5; 

s = tf('s');
w = logspace(-2, 2, 2000);

link_names = {
    '1. Интегрирующее'
    '2. Дифференцирующее'
    '3. Усилительное'
    '4. Апериодическое 1-го порядка'
    '5. Апериодическое 2-го порядка'
    '6. Колебательное'
    '7. Консервативное'
    '8. Инерционное интегрирующее'
    '9. Реальное дифференцирующее'
    '10. Форсирующее 1-го порядка'
    '11. Изодромное'
};

W = {
    K / s
    K * s
    tf(K)
    K / (T * s + 1)
    K / (T2^2 * s^2 + T1_aper * s + 1)
    K / (T2^2 * s^2 + T1 * s + 1)
    K / (T * s^2 + 1)
    K / (s * (T * s + 1))
    (K * s) / (T * s + 1)
    K * (T * s + 1)
    (K * (T * s + 1)) / s
};

figure('Name', 'АФЧХ типовых звеньев');
for i = 1:numel(W)
    subplot(4, 3, i);
    hold on; grid on;
    if i == 3
        plot(K, 0, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
        axis([0 1.5 -1 1]);
    else
        [re, im] = nyquist(W{i}, w);
        re = squeeze(re);
        im = squeeze(im);
        plot(re, im, 'b', 'LineWidth', 1.2);
        plot(re, -im, 'b--', 'LineWidth', 0.8);
    end
    switch i
        case {1, 2}, xlim([-1 1]); ylim([-6 6]);
        case 7, axis([-5 5 -1 1]);
        case 8, xlim([-2 0.5]); ylim([-6 2]);
        case 10, xlim([0 2]); ylim([-6 6]);
        case 11, xlim([0 1]); ylim([-6 6]);
    end
    title(link_names{i}, 'FontSize', 8);
    xlabel('Re', 'FontSize', 7); ylabel('Im', 'FontSize', 7);
end

figure('Name', 'ЛАЧХ типовых звеньев');
for i = 1:numel(W)
    subplot(4, 3, i);
    hold on; grid on;
    [mag, ~] = bode(W{i}, w);
    semilogx(w, 20*log10(squeeze(mag)), 'b', 'LineWidth', 1.2);
    plot([w(1) w(end)], [0 0], ':k', 'LineWidth', 0.8);
    set(gca, 'XScale', 'log');
    xlim([w(1) w(end)]);
    title(link_names{i}, 'FontSize', 8);
    xlabel('\omega, рад/с', 'FontSize', 7); ylabel('L, дБ', 'FontSize', 7);
end

figure('Name', 'ЛФЧХ типовых звеньев');
for i = 1:numel(W)
    subplot(4, 3, i);
    hold on; grid on;
    [~, phase] = bode(W{i}, w);
    semilogx(w, squeeze(phase), 'b', 'LineWidth', 1.2);
    set(gca, 'XScale', 'log');
    xlim([w(1) w(end)]);
    title(link_names{i}, 'FontSize', 8);
    xlabel('\omega, рад/с', 'FontSize', 7); ylabel('\psi, ^\circ', 'FontSize', 7);
end

K_values = [0.3, 0.5, K, 1.0, 1.4];
K_labels = compose('K = %.2g', K_values);

figure('Name', 'Влияние K на ЛАЧХ');
for i = 1:numel(W)
    subplot(4, 3, i);
    hold on; grid on;
    for k = K_values
        [mag, ~] = bode(W{i} * (k / K), w);
        semilogx(w, 20*log10(squeeze(mag)), 'LineWidth', 1.1);
    end
    plot([w(1) w(end)], [0 0], ':k', 'LineWidth', 0.8);
    set(gca, 'XScale', 'log');
    xlim([w(1) w(end)]);
    title(link_names{i}, 'FontSize', 8);
    xlabel('\omega, рад/с', 'FontSize', 7); ylabel('L, дБ', 'FontSize', 7);
    if i == 1
        legend(K_labels, 'FontSize', 6, 'Location', 'best');
    end
end

figure('Name', 'Влияние K на ЛФЧХ');
for i = 1:numel(W)
    subplot(4, 3, i);
    hold on; grid on;
    for k = K_values
        [~, phase] = bode(W{i} * (k / K), w);
        semilogx(w, squeeze(phase), 'LineWidth', 1.1);
    end
    set(gca, 'XScale', 'log');
    xlim([w(1) w(end)]);
    title(link_names{i}, 'FontSize', 8);
    xlabel('\omega, рад/с', 'FontSize', 7); ylabel('\psi, ^\circ', 'FontSize', 7);
    if i == 1
        legend(K_labels, 'FontSize', 6, 'Location', 'best');
    end
end

T_values = [0.1, 0.2, T, 0.8, 1.6];
T1_5 = [1.4, T1_aper, 2.0, 3.0, 5.0];
T2_5 = [0.2, 0.4, 0.55, T2, 0.74];
T1_6 = [0.05, 0.2, T1, 0.8, 1.2];
T2_6 = [0.35, 0.5, T2, 1.0, 1.5];

T_panels = {
    '4. Апериодическое 1-го (T)',   T_values, 'T',   @(t) K / (t*s + 1)
    '5. Апериодическое 2-го (T_1)', T1_5,     'T_1', @(t1) K / (T2^2*s^2 + t1*s + 1)
    '5. Апериодическое 2-го (T_2)', T2_5,     'T_2', @(t2) K / (t2^2*s^2 + T1_aper*s + 1)
    '6. Колебательное (T_1)',       T1_6,     'T_1', @(t1) K / (T2^2*s^2 + t1*s + 1)
    '6. Колебательное (T_2)',       T2_6,     'T_2', @(t2) K / (t2^2*s^2 + T1*s + 1)
    '7. Консервативное (T)',        T_values, 'T',   @(t) K / (t*s^2 + 1)
    '8. Инерц. интегрирующее (T)',  T_values, 'T',   @(t) K / (s*(t*s + 1))
    '9. Реальное дифф. (T)',        T_values, 'T',   @(t) (K*s) / (t*s + 1)
    '10. Форсирующее 1-го (T)',     T_values, 'T',   @(t) K * (t*s + 1)
    '11. Изодромное (T)',           T_values, 'T',   @(t) (K * (t*s + 1)) / s
};

figure('Name', 'Влияние постоянных времени на ЛАЧХ');
for j = 1:size(T_panels, 1)
    subplot(4, 3, j);
    hold on; grid on;
    vals = T_panels{j, 2};
    func = T_panels{j, 4};
    for v = vals
        [mag, ~] = bode(func(v), w);
        semilogx(w, 20*log10(squeeze(mag)), 'LineWidth', 1.1);
    end
    plot([w(1) w(end)], [0 0], ':k', 'LineWidth', 0.8);
    set(gca, 'XScale', 'log');
    xlim([w(1) w(end)]);
    title(T_panels{j, 1}, 'FontSize', 8);
    xlabel('\omega, рад/с', 'FontSize', 7); ylabel('L, дБ', 'FontSize', 7);
    legend(compose([T_panels{j, 3} ' = %.2g'], vals), 'FontSize', 6, 'Location', 'best');
end

figure('Name', 'Влияние постоянных времени на ЛФЧХ');
for j = 1:size(T_panels, 1)
    subplot(4, 3, j);
    hold on; grid on;
    vals = T_panels{j, 2};
    func = T_panels{j, 4};
    for v = vals
        [~, phase] = bode(func(v), w);
        semilogx(w, squeeze(phase), 'LineWidth', 1.1);
    end
    set(gca, 'XScale', 'log');
    xlim([w(1) w(end)]);
    title(T_panels{j, 1}, 'FontSize', 8);
    xlabel('\omega, рад/с', 'FontSize', 7); ylabel('\psi, ^\circ', 'FontSize', 7);
    legend(compose([T_panels{j, 3} ' = %.2g'], vals), 'FontSize', 6, 'Location', 'best');
end
