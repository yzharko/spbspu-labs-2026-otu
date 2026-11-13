close all; clear; clc;

pkg load control

K = 2.6;
T = 0.3;
T1 = 0.6;
T2 = 0.5;
T1_5 = 1.2;

s = tf('s');

W = {
    K / s;
    K * s;
    tf(K, 1);
    K / (T*s + 1);
    K / (T2^2*s^2 + T1_5*s + 1);
    K / (T2^2*s^2 + T1*s + 1);
    K / (T*s^2 + 1);
    K / (s*(T*s + 1));
    K*s / (T*s + 1);
    K*(T*s + 1);
    K*(T*s + 1) / s
};

omega = logspace(-2, 2, 1000);

figure(1); clf;
figure(2); clf;
figure(3); clf;

for i = 1:numel(W)
    [re, im] = nyquist(W{i}, omega);
    [mag, phase] = bode(W{i}, omega);

    re = re(:);
    im = im(:);
    mag = mag(:);
    phase = phase(:);

    L = 20 * log10(mag);

    figure(1);
    subplot(3, 4, i);

    plot(re, im, 'b', 'LineWidth', 1.5);
    hold on;
    plot(re, -im, 'b--', 'LineWidth', 1);
    hold off;

    grid on;
    title(sprintf('Звено %d', i));
    xlabel('Re');
    ylabel('Im');

    figure(2);
    subplot(3, 4, i);

    semilogx(omega, L, 'b', 'LineWidth', 1.5);

    grid on;
    title(sprintf('Звено %d', i));
    xlabel('Частота, рад/с');
    ylabel('L, дБ');

    figure(3);
    subplot(3, 4, i);

    semilogx(omega, phase, 'r', 'LineWidth', 1.5);

    grid on;
    title(sprintf('Звено %d', i));
    xlabel('Частота, рад/с');
    ylabel('Фаза, град');
end

figure(1);
set(gcf, 'Name', 'АФЧХ', 'NumberTitle', 'off');

figure(2);
set(gcf, 'Name', 'ЛАЧХ', 'NumberTitle', 'off');

figure(3);
set(gcf, 'Name', 'ЛФЧХ', 'NumberTitle', 'off');


colors = lines(5);
K_values = [K/2, K, 2*K];
K_labels = arrayfun(@(v) sprintf('K = %.2f', v), K_values, 'UniformOutput', false);

for i = 1:numel(W)
    for j = 1:numel(K_values)
        model = W{i} * K_values(j)/K;
        [mag, phase] = bode(model, omega);

        figure(4); subplot(3, 4, i); hold on;
        semilogx(omega, 20*log10(mag(:)), 'Color', colors(j,:), 'LineWidth', 1.2);
        grid on; title(sprintf('Звено %d', i));
        xlabel('Частота, рад/с'); ylabel('L, дБ');

        figure(5); subplot(3, 4, i); hold on;
        semilogx(omega, phase(:), 'Color', colors(j,:), 'LineWidth', 1.2);
        grid on; title(sprintf('Звено %d', i));
        xlabel('Частота, рад/с'); ylabel('Фаза, град');
    end
end
figure(4); subplot(3, 4, 1); legend(K_labels);
figure(5); subplot(3, 4, 1); legend(K_labels);

T_values = T * [1/5, 1/2, 1, 2, 5];
T_links = [4 7 8 9 10 11];
T_labels = arrayfun(@(v) sprintf('T = %.2f', v), T_values, 'UniformOutput', false);

for j = 1:numel(T_values)
    Ti = T_values(j);
    W_T = {
        K/(Ti*s+1),
        K/(Ti*s^2+1),
        K/(s*(Ti*s+1)),
        K*s/(Ti*s+1),
        K*(Ti*s+1),
        K*(Ti*s+1)/s
    };

    for a = 1:numel(W_T)
        [mag, phase] = bode(W_T{a}, omega);

        figure(6); subplot(2, 3, a); hold on;
        semilogx(omega, 20*log10(mag(:)), 'Color', colors(j,:), 'LineWidth', 1.3);
        grid on; title(sprintf('Звено %d', T_links(a)));
        xlabel('Частота, рад/с'); ylabel('L, дБ');

        figure(7); subplot(2, 3, a); hold on;
        semilogx(omega, phase(:), 'Color', colors(j,:), 'LineWidth', 1.3);
        grid on; title(sprintf('Звено %d', T_links(a)));
        xlabel('Частота, рад/с'); ylabel('Фаза, град');
    end
end
figure(6); subplot(2, 3, 1); legend(T_labels);
figure(7); subplot(2, 3, 1); legend(T_labels);


studies = {
    'Звено 5: T1', [1.0 1.2 1.5 1.8],       @(v) K/(T2^2*s^2 + v*s + 1);
    'Звено 5: T2', [0.30 0.40 0.45 0.50],   @(v) K/(v^2*s^2 + T1_5*s + 1);
    'Звено 6: T1', [0.18 0.30 0.60 0.90],   @(v) K/(T2^2*s^2 + v*s + 1);
    'Звено 6: T2', [0.35 0.40 0.50 0.65],   @(v) K/(v^2*s^2 + T1*s + 1)
};

for a = 1:size(studies, 1)
    title_text = studies{a, 1};
    values = studies{a, 2};
    make_model = studies{a, 3};
    labels = arrayfun(@(v) sprintf('%.2f', v), values, 'UniformOutput', false);

    for j = 1:numel(values)
        [mag, phase] = bode(make_model(values(j)), omega);

        figure(8); subplot(2, 2, a); hold on;
        semilogx(omega, 20*log10(mag(:)), 'Color', colors(j,:), 'LineWidth', 1.3);
        grid on; title(title_text);
        xlabel('Частота, рад/с'); ylabel('L, дБ');

        figure(9); subplot(2, 2, a); hold on;
        semilogx(omega, phase(:), 'Color', colors(j,:), 'LineWidth', 1.3);
        grid on; title(title_text);
        xlabel('Частота, рад/с'); ylabel('Фаза, град');
    end
    figure(8); subplot(2, 2, a); legend(labels);
    figure(9); subplot(2, 2, a); legend(labels);
end
