clc;
clear;
close all;

K = 2.1;
T = 0.9;
T1 = 1.0;
T2 = 1.1;
T1_ap = 3.0;

s = tf('s');
w = logspace(-2, 2, 1500);
names = {'Интегрирующее', 'Идеальное дифференцирующее', ...
    'Усилительное', 'Апериодическое 1-го порядка', ...
    'Апериодическое 2-го порядка', 'Колебательное', ...
    'Консервативное', 'Инерционное интегрирующее', ...
    'Реальное дифференцирующее', 'Форсирующее 1-го порядка', ...
    'Изодромное'};
W = make_links(s, K, T, T1, T2, T1_ap);

out = fullfile(fileparts(mfilename('fullpath')), 'L3_graphs', 'combined');
if ~exist(out, 'dir')
    mkdir(out);
end

titles = cell(1, 11);
for n = 1:11
    titles{n} = sprintf('%d. %s', n, names{n});
end
show_bode(W(:), w, titles, {}, out, 'base', 'Все звенья');

figure('Name', 'АФЧХ всех звеньев', 'NumberTitle', 'off', ...
    'Position', [50 50 1500 1100]);
tiledlayout(4, 3, 'TileSpacing', 'compact');
for n = 1:11
    nexttile;
    [re, im] = nyquist(W{n}, w);
    re = squeeze(re);
    im = squeeze(im);
    if max(re) - min(re) < 1e-10
        re(:) = round(mean(re), 10);
    end
    im(abs(im) < 1e-10) = 0;
    if n == 3
        plot(K, 0, 'bo', 'MarkerFaceColor', 'b');
    elseif n == 7
        re([false; diff(sign(re)) ~= 0]) = NaN;
        plot(re, im, 'b', 'LineWidth', 1.1);
        xlim([-15 15]);
        ylim([-1 1]);
    else
        plot(re, im, 'b', re, -im, 'b', 'LineWidth', 1.1);
    end
    if max(re) - min(re) < 1e-10
        xlim([re(1)-1 re(1)+1]);
    end
    grid on;
    xlabel('Re W');
    ylabel('Im W');
    title(titles{n}, 'FontSize', 9);
end
sgtitle('АФЧХ всех звеньев');
exportgraphics(gcf, fullfile(out, 'nyquist_all.png'));

K_values = K*[0.5 1 2];
modelsK = cell(11, 3);
labelsK = repmat({compose('K = %g', K_values)}, 11, 1);
for n = 1:11
    for j = 1:3
        modelsK{n, j} = W{n}*K_values(j)/K;
    end
end
show_bode(modelsK, w, titles, labelsK, out, 'K', 'Влияние K');

T_values = T*[0.5 1 2];
numbers = [4 7 8 9 10 11];
modelsT = cell(10, 3);
labelsT = cell(10, 1);
titlesT = cell(10, 1);
for i = 1:6
    n = numbers(i);
    for j = 1:3
        V = make_links(s, K, T_values(j), T1, T2, T1_ap);
        modelsT{i, j} = V{n};
    end
    labelsT{i} = compose('T = %g', T_values);
    titlesT{i} = titles{n};
end

row = 6;
for n = [5 6]
    if n == 5
        T1_values = [2.4 3.0 4.5];
        T2_values = [0.7 1.1 1.4];
        baseT1 = T1_ap;
    else
        T1_values = [0.5 1.0 1.5];
        T2_values = [0.7 1.1 1.5];
        baseT1 = T1;
    end
    for parameter = 1:2
        row = row + 1;
        for j = 1:3
            if parameter == 1
                modelsT{row, j} = K/(T2^2*s^2 + T1_values(j)*s + 1);
            else
                modelsT{row, j} = K/(T2_values(j)^2*s^2 + baseT1*s + 1);
            end
        end
        if parameter == 1
            labelsT{row} = compose('T1 = %g', T1_values);
        else
            labelsT{row} = compose('T2 = %g', T2_values);
        end
        titlesT{row} = sprintf('%d. Влияние T%d', n, parameter);
    end
end
show_bode(modelsT, w, titlesT, labelsT, out, 'T', 'Влияние постоянных времени');

fprintf('Низкочастотное усиление: 20*lg(K) = %.3f дБ\n', 20*log10(K));
fprintf('При удвоении K сдвиг ЛАЧХ равен %.3f дБ\n', 20*log10(2));
for j = 1:3
    fprintf('T = %.2f: 1/T = %.4f рад/с; 1/sqrt(T) = %.4f рад/с\n', ...
        T_values(j), 1/T_values(j), 1/sqrt(T_values(j)));
end
tau = [(T1_ap + sqrt(T1_ap^2 - 4*T2^2))/2, ...
       (T1_ap - sqrt(T1_ap^2 - 4*T2^2))/2];
fprintf('Звено 5: сопрягающие частоты %.4f и %.4f рад/с\n', 1/tau(1), 1/tau(2));
zeta = T1/(2*T2);
fprintf('Звено 6: w0 = %.4f рад/с; zeta = %.4f\n', 1/T2, zeta);
if zeta < 1/sqrt(2)
    wr = sqrt(1 - 2*zeta^2)/T2;
    Ar = K/(2*zeta*sqrt(1 - zeta^2));
    fprintf('Резонанс: wr = %.4f рад/с; Amax = %.4f; Lmax = %.3f дБ\n', ...
        wr, Ar, 20*log10(Ar));
end
fprintf('Графики сохранены в %s\n', out);

function W = make_links(s, K, T, T1, T2, T1_ap)
W = {K/s, K*s, tf(K), K/(T*s + 1), ...
    K/(T2^2*s^2 + T1_ap*s + 1), ...
    K/(T2^2*s^2 + T1*s + 1), ...
    K/(T*s^2 + 1), K/(s*(T*s + 1)), ...
    K*s/(T*s + 1), K*(T*s + 1), K*(T*s + 1)/s};
end

function show_bode(models, w, titles, labels, out, prefix, name)
ylabels = {'L, дБ', '\phi, град'};
chartNames = {'ЛАЧХ', 'ЛФЧХ'};
for mode = 1:2
    figure('Name', [name ' — ' chartNames{mode}], 'NumberTitle', 'off', ...
        'Position', [50 50 1500 1100]);
    tiledlayout(4, 3, 'TileSpacing', 'compact');
    for n = 1:size(models, 1)
        ax = nexttile;
        hold on;
        for j = 1:size(models, 2)
            [mag, phase] = bode(models{n, j}, w);
            if mode == 1
                semilogx(w, 20*log10(squeeze(mag)), 'LineWidth', 1.1);
            else
                semilogx(w, squeeze(phase), 'LineWidth', 1.1);
            end
        end
        set(gca, 'XScale', 'log');
        grid on;
        xlabel('\omega, рад/с');
        ylabel(ylabels{mode});
        title(titles{n}, 'FontSize', 9);
        if ~isempty(labels)
            legend(ax, labels{n}, 'Location', 'northeast', 'FontSize', 8);
        end
    end
    sgtitle([name ' — ' chartNames{mode}]);
    drawnow;
    exportgraphics(gcf, fullfile(out, sprintf('%s_%d.png', prefix, mode)));
end
end
