clear; clc; close all;

K = 0.9;
T = 0.6;
T1 = 0.8;
T2 = 0.7;
T1a = 1.6;

w = logspace(-2, 2, 1500);
names = {'Интегрирующее', 'Идеальное дифференцирующее', ...
    'Усилительное', 'Апериодическое 1-го порядка', ...
    'Апериодическое 2-го порядка', 'Колебательное', ...
    'Консервативное', 'Инерционное интегрирующее', ...
    'Реальное дифференцирующее', 'Форсирующее 1-го порядка', ...
    'Изодромное'};

W = make_links(K, T, T1, T2, T1a);

show_bode(W, names, w, 1:11, {}, 'Исходные параметры');

for first = 1:3:11
    indices = first:min(first + 2, 11);
    figure('Name', sprintf('АФЧХ, звенья %s', num2str(indices)), ...
        'Color', 'w', 'Position', [60 100 1400 500]);
    for column = 1:numel(indices)
        n = indices(column);
        [re, im] = nyquist(W{n}, w);
        re = squeeze(re);
        im = squeeze(im);
        switch n
            case {1, 2}
                re(:) = 0;
            case 7
                im(:) = 0;
                cut = find(w >= 1/sqrt(T), 1);
                re(cut) = NaN;
                im(cut) = NaN;
            case 10
                re(:) = K;
            case 11
                re(:) = K*T;
        end

        subplot(1, numel(indices), column);
        plot(re, im, 'b', re, -im, 'r--', 'LineWidth', 1.5);
        hold on; grid on; box on;
        if n == 3
            plot(K, 0, 'bo', 'MarkerFaceColor', 'b');
        elseif n == 7
            xlim([-5 5]); ylim([-1 1]);
        end
        title(sprintf('%d. %s', n, names{n}), 'FontSize', 12, 'Color', 'k');
        xlabel('Re W(j\omega)', 'Color', 'k');
        ylabel('Im W(j\omega)', 'Color', 'k');
        set(gca, 'FontSize', 12, 'Color', 'w', 'XColor', 'k', ...
            'YColor', 'k', 'Toolbar', []);
        legend('\omega > 0', '\omega < 0', 'Location', 'best', ...
            'TextColor', 'k', 'Color', 'w');
    end
    sgtitle('Амплитудно-фазовые характеристики', 'Color', 'k');
end

K_values = K*[0.5 1 2];
WK = cell(numel(K_values), 11);
for j = 1:numel(K_values)
    WK(j, :) = make_links(K_values(j), T, T1, T2, T1a);
end
show_bode(WK, names, w, 1:11, ...
    compose('K = %.2g', K_values), 'Влияние K');

fprintf('Влияние K на уровень ЛАЧХ\n');
fprintf('K       Смещение относительно K = %.2f, дБ\n', K);
fprintf('%.2f    %+.3f\n', [K_values; 20*log10(K_values/K)]);

T_values = T*[0.5 1 2];
WT = cell(numel(T_values), 11);
for j = 1:numel(T_values)
    WT(j, :) = make_links(K, T_values(j), T1, T2, T1a);
end
show_bode(WT, names, w, [4 7 8 9 10 11], ...
    compose('T = %.2g', T_values), 'Влияние T');

fprintf('\nХарактерные частоты при изменении T\n');
fprintf('T       1/T, рад/с    1/sqrt(T), рад/с\n');
fprintf('%.2f    %.4f        %.4f\n', ...
    [T_values; 1./T_values; 1./sqrt(T_values)]);

T1a_values = [1.4 1.6 2.0];
T1_values = [0.4 0.8 1.2];
T2_values = [0.5 0.7 0.8];
WT1 = cell(3, 11);
WT2 = cell(3, 11);
for j = 1:3
    WT1(j, :) = make_links(K, T, T1_values(j), T2, T1a_values(j));
    WT2(j, :) = make_links(K, T, T1, T2_values(j), T1a);
end
show_bode(WT1, names, w, 5, compose('T_1 = %.2g', T1a_values), ...
    'Влияние T1 при T2 = 0.7');
show_bode(WT1, names, w, 6, compose('T_1 = %.2g', T1_values), ...
    'Влияние T1 при T2 = 0.7');
show_bode(WT2, names, w, [5 6], compose('T_2 = %.2g', T2_values), ...
    'Влияние T2 при постоянном T1');

tau = [(T1a + sqrt(T1a^2 - 4*T2^2))/2, ...
       (T1a - sqrt(T1a^2 - 4*T2^2))/2];
zeta = T1/(2*T2);
fprintf('\nИсходные параметры звеньев второго порядка\n');
fprintf('Звено 5: сопрягающие частоты %.4f и %.4f рад/с\n', 1./tau);
fprintf('Звено 6: omega0 = %.4f рад/с, zeta = %.4f\n', 1/T2, zeta);

function W = make_links(K, T, T1, T2, T1a)
    s = tf('s');
    W = {K/s, K*s, tf(K), K/(T*s + 1), ...
        K/(T2^2*s^2 + T1a*s + 1), ...
        K/(T2^2*s^2 + T1*s + 1), K/(T*s^2 + 1), ...
        K/(s*(T*s + 1)), K*s/(T*s + 1), ...
        K*(T*s + 1), K*(T*s + 1)/s};
end

function show_bode(models, names, w, indices, labels, caption)
    for first = 1:3:numel(indices)
        group = indices(first:min(first + 2, numel(indices)));
        columns = numel(group);
        figure('Name', sprintf('%s, звенья %s', caption, num2str(group)), ...
            'Color', 'w', 'Position', [60 60 1400 800]);
        for column = 1:columns
            n = group(column);
            ax_mag = subplot(2, columns, column);
            hold on; grid on; box on;
            ax_phase = subplot(2, columns, columns + column);
            hold on; grid on; box on;

            for j = 1:size(models, 1)
                [mag, phase] = bode(models{j, n}, w);
                mag = squeeze(mag);
                phase = mod(squeeze(phase) + 180, 360) - 180;
                if n == 7
                    [~, den] = tfdata(models{j, n}, 'v');
                    cut = find(w >= sqrt(den(end)/den(1)), 1);
                    mag(cut) = NaN;
                    phase(cut) = NaN;
                end
                semilogx(ax_mag, w, 20*log10(mag), 'LineWidth', 1.5);
                semilogx(ax_phase, w, phase, 'LineWidth', 1.5);
            end

            title(ax_mag, sprintf('%d. %s', n, names{n}), ...
                'FontSize', 12, 'Color', 'k');
            xlabel(ax_mag, '\omega, рад/с', 'Color', 'k');
            ylabel(ax_mag, 'L(\omega), дБ', 'Color', 'k');
            xlabel(ax_phase, '\omega, рад/с', 'Color', 'k');
            ylabel(ax_phase, '\psi(\omega), град', 'Color', 'k');
            yline(ax_mag, 0, ':k', 'HandleVisibility', 'off');
            set([ax_mag ax_phase], 'FontSize', 12, 'XScale', 'log', ...
                'XLim', [w(1) w(end)], 'Color', 'w', 'XColor', 'k', ...
                'YColor', 'k', 'Toolbar', []);
            if ~isempty(labels)
                legend(ax_mag, labels, 'Location', 'best', ...
                    'TextColor', 'k', 'Color', 'w');
            end
        end
        sgtitle(caption, 'Color', 'k');
    end
end
