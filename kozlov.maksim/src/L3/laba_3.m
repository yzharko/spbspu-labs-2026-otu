clc;
clear;
close all;

s = tf('s');
K_base   = 4.0;   
T_base   = 1.5;   
T1_base  = 1.6;   
T2_base  = 2.0;   
T1a_base = 4.5;   


names = {'1. Интегрирующее', '2. Дифференцирующее', '3. Усилительное', ...
         '4. Апериодическое 1-го', '5. Апериодическое 2-го', '6. Колебательное', ...
         '7. Консервативное', '8. Инерционное интегр.', '9. Реальное дифф.', ...
         '10. Форсирующее', '11. Изодромное'};

elements = @(K, T, T1, T2, T1a) {K/s, ...
                                 K*s, ...
                                 tf(K), ...
                                 K/(T*s + 1), ...
                                 K/(T2^2*s^2 + T1a*s + 1), ...
                                 K/(T2^2*s^2 + T1*s + 1), ...
                                 K/(T*s^2 + 1), ...
                                 K/(s*(T*s + 1)), ...
                                 K*s/(T*s + 1), ...
                                 K*(T*s + 1), ...
                                 K*(T*s + 1)/s};

W = elements(K_base, T_base, T1_base, T2_base, T1a_base);

folder = fileparts(mfilename('fullpath'));
if isempty(folder)
    folder = pwd;
end

w  = logspace(-2, 2, 2000);
wn = logspace(-1.3, 1.3, 4000);
w2 = logspace(-2, 2, 1000);

f1 = figure('Color', 'w', 'Position', [50 50 1400 900]);
for i = 1:numel(names)
    [re, im] = nyquist(W{i}, wn);
    re = squeeze(re);
    im = squeeze(im);
    if isscalar(re)
        re = repmat(re, size(wn));
        im = repmat(im, size(wn));
    end

    subplot(3,4,i); hold on; grid on;
    plot(re, im, 'b', 'LineWidth', 2);
    plot(re(1),   im(1),   'go', 'MarkerFaceColor', 'g');
    plot(re(end), im(end), 'ks', 'MarkerFaceColor', 'k');

    if max(re) - min(re) < 1e-3, xlim(round(mean(re)) + [-1 1]); end
    if i == 3, xlim([2 6]);   ylim([-1 1]); end
    if i == 7, xlim([-5 15]); ylim([-1 1]); end
    if i == 1, legend({'АФЧХ', '\omega=0.05', '\omega=20'}, 'Location', 'best'); end

    title(names{i}, 'FontWeight', 'bold');
    xlabel('Re');
    ylabel('Im');
end
exportgraphics(f1, fullfile(folder, '01_АФЧХ.png'));

f2 = figure('Color', 'w', 'Position', [80 80 1600 1000]);
f3 = figure('Color', 'w', 'Position', [110 110 1600 1000]);
for i = 1:numel(names)
    [mag, ph] = bode(W{i}, w);

    figure(f2); subplot(3,4,i); hold on; grid on;
    semilogx(w, 20*log10(squeeze(mag)), 'b', 'LineWidth', 2);
    set(gca, 'XScale', 'log');
    title(names{i}, 'FontWeight', 'bold');
    xlabel('\omega, рад/с');
    ylabel('L, дБ');

    figure(f3); subplot(3,4,i); hold on; grid on;
    semilogx(w, squeeze(ph), 'm', 'LineWidth', 2);
    set(gca, 'XScale', 'log');
    title(names{i}, 'FontWeight', 'bold');
    xlabel('\omega, рад/с');
    ylabel('\psi, °');
end
exportgraphics(f2, fullfile(folder, '02_ЛАЧХ.png'));
exportgraphics(f3, fullfile(folder, '03_ЛФЧХ.png'));

K_values = [1 2 4 6 8];
T_values = [0.2 0.5 1 2 5];
styles = {'-', '--', ':', '-.', '-'};
widths = [3 2.5 2 2 1.5];

WK = cell(1, length(K_values));
labelsK = cell(1, length(K_values));
for j = 1:length(K_values)
    WK{j} = elements(K_values(j), T_base, T1_base, T2_base, T1a_base);
    labelsK{j} = sprintf('K = %.2g', K_values(j));
end

f4 = figure('Color', 'w', 'Position', [30 30 1800 1000]);
f5 = figure('Color', 'w', 'Position', [30 30 1800 1000]);
for i = 1:numel(names)
    for j = 1:length(K_values)
        [mag, ph] = bode(WK{j}{i}, w2);

        figure(f4); subplot(3,4,i); hold on; grid on;
        semilogx(w2, 20*log10(squeeze(mag)), 'LineStyle', styles{j}, 'LineWidth', widths(j));

        figure(f5); subplot(3,4,i); hold on; grid on;
        semilogx(w2, squeeze(ph), 'LineStyle', styles{j}, 'LineWidth', widths(j));
    end

    figure(f4); subplot(3,4,i);
    set(gca, 'XScale', 'log');
    title(names{i}, 'FontWeight', 'bold');
    xlabel('\omega, рад/с');
    ylabel('L, дБ');
    if i == 1, legend(labelsK, 'Location', 'best'); end

    figure(f5); subplot(3,4,i);
    set(gca, 'XScale', 'log');
    yl = ylim;
    if diff(yl) < 20, ylim(mean(yl) + [-90 90]); end
    title(names{i}, 'FontWeight', 'bold');
    xlabel('\omega, рад/с');
    ylabel('\psi, °');
    if i == 1, legend(labelsK, 'Location', 'best'); end
end
figure(f4); sgtitle('Влияние K на ЛАЧХ', 'FontWeight', 'bold');
figure(f5); sgtitle('Влияние K на ЛФЧХ', 'FontWeight', 'bold');
exportgraphics(f4, fullfile(folder, '04_Влияние_K_ЛАЧХ.png'));
exportgraphics(f5, fullfile(folder, '05_Влияние_K_ЛФЧХ.png'));

WT = cell(1, length(T_values));
labelsT = cell(1, length(T_values));
for j = 1:length(T_values)
    m = T_values(j);
    Wm = elements(K_base, T_base*m, T1_base*m, T2_base*m, T1a_base*m);
    w0 = 1/sqrt(T_base*m);
    W7 = K_base/(s^2/w0^2 + 1);
    WT{j} = [Wm(1:6), {W7}, Wm(8:11)];
    labelsT{j} = sprintf('T x %.2g', m);
end

f6 = figure('Color', 'w', 'Position', [30 30 1800 1000]);
f7 = figure('Color', 'w', 'Position', [30 30 1800 1000]);
for i = 1:numel(names)
    for j = 1:length(T_values)
        [mag, ph] = bode(WT{j}{i}, w2);

        figure(f6); subplot(3,4,i); hold on; grid on;
        semilogx(w2, 20*log10(squeeze(mag)), 'LineStyle', styles{j}, 'LineWidth', widths(j));

        figure(f7); subplot(3,4,i); hold on; grid on;
        semilogx(w2, squeeze(ph), 'LineStyle', styles{j}, 'LineWidth', widths(j));
    end

    figure(f6); subplot(3,4,i);
    set(gca, 'XScale', 'log');
    title(names{i}, 'FontWeight', 'bold');
    xlabel('\omega, рад/с');
    ylabel('L, дБ');
    if i == 1, legend(labelsT, 'Location', 'best'); end

    figure(f7); subplot(3,4,i);
    set(gca, 'XScale', 'log');
    yl = ylim;
    if diff(yl) < 20, ylim(mean(yl) + [-90 90]); end
    title(names{i}, 'FontWeight', 'bold');
    xlabel('\omega, рад/с');
    ylabel('\psi, °');
    if i == 1, legend(labelsT, 'Location', 'best'); end
end
figure(f6); sgtitle('Влияние T на ЛАЧХ', 'FontWeight', 'bold');
figure(f7); sgtitle('Влияние T на ЛФЧХ', 'FontWeight', 'bold');
exportgraphics(f6, fullfile(folder, '06_Влияние_T_ЛАЧХ.png'));
exportgraphics(f7, fullfile(folder, '07_Влияние_T_ЛФЧХ.png'));

