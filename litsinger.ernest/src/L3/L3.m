clc;
clear;
close all;

s = tf('s');

K   = 0.8;
T   = 0.5;
T1  = 0.7;
T2  = 0.5;
T1a = 1.5;

names = {
    '1. Интегрирующее';
    '2. Дифференцирующее';
    '3. Усилительное';
    '4. Апериодическое 1-го пор.';
    '5. Апериодическое 2-го пор.';
    '6. Колебательное';
    '7. Консервативное';
    '8. Инерционное интегр.';
    '9. Реальное дифференц.';
    '10. Форсирующее';
    '11. Изодромное'};

links = {
    @(K, T, T1, T2, T1a) K/s;
    @(K, T, T1, T2, T1a) K*s;
    @(K, T, T1, T2, T1a) tf(K);
    @(K, T, T1, T2, T1a) K/(T*s + 1);
    @(K, T, T1, T2, T1a) K/(T2^2*s^2 + T1a*s + 1);
    @(K, T, T1, T2, T1a) K/(T2^2*s^2 + T1*s + 1);
    @(K, T, T1, T2, T1a) K/(T*s^2 + 1);
    @(K, T, T1, T2, T1a) K/(s*(T*s + 1));
    @(K, T, T1, T2, T1a) K*s/(T*s + 1);
    @(K, T, T1, T2, T1a) K*(T*s + 1);
    @(K, T, T1, T2, T1a) K*(T*s + 1)/s};

idxAp2  = 5;
idxOsc  = 6;
idxCons = 7;

nLinks = numel(links);
nCols  = 4;
nRows  = ceil(nLinks / nCols);
lw     = 1.2;

w  = logspace(-2, 3, 2000);
wa = logspace(-1, 2, 2000);
consAxis = [-5 5 -1 1];

K_values  = K  * [0.25 0.5 1 2 4];
T_values  = T  * [0.25 0.5 1 2 4];
T2_ap2    = [0.2 0.3 0.5 0.6 0.75];
T2_osc    = [0.4 0.5 0.7 1 1.5];
T1_osc    = [0.1 0.3 0.5 0.7 0.9];
nVar      = numel(K_values);

figure('Name', 'АФЧХ');
for n = 1:nLinks
    W = links{n}(K, T, T1, T2, T1a);
    [re, im] = nyquist(W, wa);
    re = squeeze(re);
    im = squeeze(im);
    subplot(nRows, nCols, n);
    plot(re, im, 'LineWidth', lw);
    hold on;
    plot(re(1), im(1), 'ko', 'MarkerFaceColor', 'k');
    grid on;
    title(names{n});
    xlabel('U(\omega)');
    ylabel('V(\omega)');
    if n == idxCons
        axis(consAxis);
    end
end

figure('Name', 'ЛАЧХ');
for n = 1:nLinks
    [L, ~] = lfc(links{n}(K, T, T1, T2, T1a), w);
    subplot(nRows, nCols, n);
    semilogx(w, L, 'LineWidth', lw);
    grid on;
    title(names{n});
    xlabel('\omega, рад/с');
    ylabel('L, дБ');
end

figure('Name', 'ЛФЧХ');
for n = 1:nLinks
    [~, phi] = lfc(links{n}(K, T, T1, T2, T1a), w);
    subplot(nRows, nCols, n);
    semilogx(w, phi, 'LineWidth', lw);
    grid on;
    title(names{n});
    xlabel('\omega, рад/с');
    ylabel('\psi, град');
end

figure('Name', 'Влияние K на ЛАЧХ');
for n = 1:nLinks
    subplot(nRows, nCols, n);
    hold on;
    for i = 1:nVar
        [L, ~] = lfc(links{n}(K_values(i), T, T1, T2, T1a), w);
        semilogx(w, L, 'LineWidth', lw);
    end
    set(gca, 'XScale', 'log');
    grid on;
    title(names{n});
    xlabel('\omega, рад/с');
    ylabel('L, дБ');
end
legend(compose('K = %.2g', K_values), 'Location', 'best');

figure('Name', 'Влияние K на ЛФЧХ');
for n = 1:nLinks
    subplot(nRows, nCols, n);
    hold on;
    for i = 1:nVar
        [~, phi] = lfc(links{n}(K_values(i), T, T1, T2, T1a), w);
        semilogx(w, phi, 'LineWidth', lw);
    end
    set(gca, 'XScale', 'log');
    grid on;
    title(names{n});
    xlabel('\omega, рад/с');
    ylabel('\psi, град');
end
legend(compose('K = %.2g', K_values), 'Location', 'best');

WT = cell(nLinks, nVar);
for n = 1:nLinks
    for i = 1:nVar
        if n == idxAp2
            WT{n, i} = links{n}(K, T, T1, T2_ap2(i), T1a);
        elseif n == idxOsc
            WT{n, i} = links{n}(K, T, T1, T2_osc(i), T1a);
        else
            WT{n, i} = links{n}(K, T_values(i), T1, T2, T1a);
        end
    end
end

figure('Name', 'Влияние T на ЛАЧХ');
for n = 1:nLinks
    subplot(nRows, nCols, n);
    hold on;
    for i = 1:nVar
        [L, ~] = lfc(WT{n, i}, w);
        semilogx(w, L, 'LineWidth', lw);
    end
    set(gca, 'XScale', 'log');
    grid on;
    title(names{n});
    xlabel('\omega, рад/с');
    ylabel('L, дБ');
    if n == idxAp2
        legend(compose('T2 = %.2g', T2_ap2), 'Location', 'best');
    elseif n == idxOsc
        legend(compose('T2 = %.2g', T2_osc), 'Location', 'best');
    end
end
legend(compose('T = %.3g', T_values), 'Location', 'best');

figure('Name', 'Влияние T на ЛФЧХ');
for n = 1:nLinks
    subplot(nRows, nCols, n);
    hold on;
    for i = 1:nVar
        [~, phi] = lfc(WT{n, i}, w);
        semilogx(w, phi, 'LineWidth', lw);
    end
    set(gca, 'XScale', 'log');
    grid on;
    title(names{n});
    xlabel('\omega, рад/с');
    ylabel('\psi, град');
end

figure('Name', 'Колебательное звено, влияние T1');
for i = 1:numel(T1_osc)
    [L, phi] = lfc(links{idxOsc}(K, T, T1_osc(i), T2, T1a), w);
    subplot(1, 2, 1);
    semilogx(w, L, 'LineWidth', lw);
    hold on;
    subplot(1, 2, 2);
    semilogx(w, phi, 'LineWidth', lw);
    hold on;
end
subplot(1, 2, 1);
grid on;
title('ЛАЧХ');
xlabel('\omega, рад/с');
ylabel('L, дБ');
legend(compose('T1 = %.1f', T1_osc), 'Location', 'best');
subplot(1, 2, 2);
grid on;
title('ЛФЧХ');
xlabel('\omega, рад/с');
ylabel('\psi, град');

fprintf('Сопрягающая частота 1/T = %.2f рад/с\n', 1/T);
fprintf('Собственная частота колебательного звена 1/T2 = %.2f рад/с\n', 1/T2);
fprintf('Коэффициент затухания T1/(2*T2) = %.2f\n', T1/(2*T2));
fprintf('Резонансная частота консервативного звена 1/sqrt(T) = %.3f рад/с\n', 1/sqrt(T));
fprintf('20*lg(K) = %.2f дБ\n', 20*log10(K));

function [L, phi] = lfc(W, w)
    [mag, phase] = bode(W, w);
    L   = 20*log10(squeeze(mag));
    phi = squeeze(phase);
    if isscalar(L)
        L   = L * ones(size(w));
        phi = phi * ones(size(w));
    end
end
