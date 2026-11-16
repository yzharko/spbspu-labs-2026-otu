clear; clc; close all;

K0 = 3;
T0 = 1;
T1 = 1.2;
T2 = 1.5;
T1a = 3.75;

syms s
syms w K T positive

names = {'Интегрирующее', 'Дифференцирующее', 'Усилительное', ...
    'Апериодическое 1-го порядка', 'Апериодическое 2-го порядка', 'Колебательное', ...
    'Консервативное', 'Инерционное интегрирующее', 'Реальное дифференцирующее', ...
    'Форсирующее 1-го порядка', 'Изодромное'};
W = [K/s; K*s; K; K/(T*s + 1); K/(T2^2*s^2 + T1a*s + 1); K/(T2^2*s^2 + T1*s + 1); ...
    K/(T*s^2 + 1); K/(s*(T*s + 1)); K*s/(T*s + 1); K*(T*s + 1); K*(T*s + 1)/s];
N = numel(names);

W_jw = subs(W, s, 1i*w);
U = simplify(real(W_jw));
V = simplify(imag(W_jw));
A = simplify(sqrt(U.^2 + V.^2));
for i = 1:N
    fprintf('%d. %s\n   U = %s\n   V = %s\n   A = %s\n', i, names{i}, ...
        char(U(i)), char(V(i)), char(A(i)));
end

s_tf = tf('s');
W_tf = {K0/s_tf, K0*s_tf, tf(K0), K0/(T0*s_tf + 1), K0/(T2^2*s_tf^2 + T1a*s_tf + 1), ...
    K0/(T2^2*s_tf^2 + T1*s_tf + 1), K0/(T0*s_tf^2 + 1), K0/(s_tf*(T0*s_tf + 1)), ...
    K0*s_tf/(T0*s_tf + 1), K0*(T0*s_tf + 1), K0*(T0*s_tf + 1)/s_tf};

omega = logspace(-2, 2, 1000);
afc_limit = 10;
afc_margin = 1;
W_jw_value = subs(W_jw, [K T], [K0 T0]);
dev = zeros(N, 2);

fig_names = {'АФЧХ', 'ЛАЧХ', 'ЛФЧХ'};
x_labels = {'U(\omega)', '\omega, рад/с', '\omega, рад/с'};
y_labels = {'V(\omega)', 'L(\omega), дБ', '\psi(\omega), град'};
for k = 1:3
    figure('Name', fig_names{k}, 'Position', [50 50 1000 1200]);
    layouts(k) = tiledlayout(4, 3, 'TileSpacing', 'compact');
    title(layouts(k), sprintf('%s типовых звеньев при K = %g, T = %g с', fig_names{k}, K0, T0));
end
for i = 1:N
    [mag, phase] = bode(W_tf{i}, omega);
    mag = squeeze(mag);
    phase = squeeze(phase);
    [re, im] = nyquist(W_tf{i}, omega);
    re = squeeze(re);
    im = squeeze(im);
    W_value = double(subs(W_jw_value(i), w, omega)).';
    dev(i, :) = [max(abs(re + 1i*im - W_value)./abs(W_value)), ...
        max(abs(mag.*exp(1i*deg2rad(phase)) - W_value)./abs(W_value))];
    re(mag > afc_limit) = NaN;
    im(mag > afc_limit) = NaN;

    axes_i = [nexttile(layouts(1), i), nexttile(layouts(2), i), nexttile(layouts(3), i)];
    plot(axes_i(1), re, im, 'LineWidth', 1.5);
    axis(axes_i(1), [min([re; 0]) max([re; 0]) min([im; 0]) max([im; 0])] + ...
        afc_margin*[-1 1 -1 1]);
    if isstatic(W_tf{i})
        hold(axes_i(1), 'on');
        plot(axes_i(1), re(1), im(1), 'o', 'LineWidth', 2);
    end
    semilogx(axes_i(2), omega, 20*log10(mag), 'LineWidth', 1.5);
    semilogx(axes_i(3), omega, phase, 'LineWidth', 1.5);
    for k = 1:3
        grid(axes_i(k), 'on');
        title(axes_i(k), sprintf('%d. %s', i, names{i}));
        xlabel(axes_i(k), x_labels{k});
        ylabel(axes_i(k), y_labels{k});
    end
end
fprintf('\nНаибольшее относительное расхождение аналитики с nyquist: %.1e, с bode: %.1e\n', ...
    max(dev));

K_values = K0*[0.5 1 2];
crossing = NaN(N, numel(K_values));
phase_change = 0;
figure('Name', 'Влияние K', 'Position', [50 50 1000 1200]);
layout_K = tiledlayout(4, 3, 'TileSpacing', 'compact');
title(layout_K, sprintf('ЛАЧХ типовых звеньев при T = %g с и различных K', T0));
for i = 1:N
    ax = nexttile(layout_K, i);
    hold(ax, 'on');
    [~, phase_base] = bode(W_tf{i}, omega);
    for j = 1:numel(K_values)
        [mag, phase] = bode(K_values(j)/K0*W_tf{i}, omega);
        L = 20*log10(squeeze(mag));
        plot(ax, omega, L, 'LineWidth', 1.5, 'DisplayName', sprintf('K = %g', K_values(j)));
        idx = find(diff(sign(L)) ~= 0, 1);
        if ~isempty(idx)
            crossing(i, j) = interp1(L(idx:idx + 1), omega(idx:idx + 1), 0);
        end
        phase_change = max(phase_change, max(abs(phase - phase_base), [], 'all'));
    end
    set(ax, 'XScale', 'log');
    grid(ax, 'on');
    title(ax, sprintf('%d. %s', i, names{i}));
    xlabel(ax, '\omega, рад/с');
    ylabel(ax, 'L(\omega), дБ');
end
legend_K = legend(ax);
legend_K.Layout.Tile = N + 1;
fprintf('\nЧастота пересечения ЛАЧХ с осью 0 дБ, рад/с, при K = %g; %g; %g\n', K_values);
for i = 1:N
    fprintf('%2d %8.3f %8.3f %8.3f\n', i, crossing(i, :));
end
fprintf('Наибольшее изменение ЛФЧХ при изменении K: %.1e град\n', phase_change);

T_scale = [0.5 1 2];
T_study = { ...
    4, 'T', T0*T_scale, @(value) K0/(value*s_tf + 1); ...
    5, 'T_1', [2*T2 T1a 2*T1a], @(value) K0/(T2^2*s_tf^2 + value*s_tf + 1); ...
    5, 'T_2', [T2/2 T2 T1a/2], @(value) K0/(value^2*s_tf^2 + T1a*s_tf + 1); ...
    6, 'T_1', T1*T_scale, @(value) K0/(T2^2*s_tf^2 + value*s_tf + 1); ...
    6, 'T_2', T2*T_scale, @(value) K0/(value^2*s_tf^2 + T1*s_tf + 1); ...
    7, 'T', T0*T_scale, @(value) K0/(value*s_tf^2 + 1); ...
    8, 'T', T0*T_scale, @(value) K0/(s_tf*(value*s_tf + 1)); ...
    9, 'T', T0*T_scale, @(value) K0*s_tf/(value*s_tf + 1); ...
    10, 'T', T0*T_scale, @(value) K0*(value*s_tf + 1); ...
    11, 'T', T0*T_scale, @(value) K0*(value*s_tf + 1)/s_tf};
fig_T_names = {'ЛАЧХ', 'ЛФЧХ'};
for k = 1:2
    figure('Name', ['Влияние T: ' fig_T_names{k}], 'Position', [50 50 1000 1200]);
    layouts_T(k) = tiledlayout(4, 3, 'TileSpacing', 'compact');
    title(layouts_T(k), sprintf('%s типовых звеньев при K = %g и различных постоянных времени, с', ...
        fig_T_names{k}, K0));
end
fprintf('\nСопрягающие частоты, рад/с\n');
for k = 1:size(T_study, 1)
    [link, param, values, make_W] = T_study{k, :};
    axes_k = [nexttile(layouts_T(1), k), nexttile(layouts_T(2), k)];
    hold(axes_k(1), 'on');
    hold(axes_k(2), 'on');
    for value = values
        W_T = make_W(value);
        [mag, phase] = bode(W_T, omega);
        label = sprintf('%s = %g', param, value);
        plot(axes_k(1), omega, 20*log10(squeeze(mag)), 'LineWidth', 1.5, 'DisplayName', label);
        plot(axes_k(2), omega, squeeze(phase), 'LineWidth', 1.5, 'DisplayName', label);
        roots_abs = abs([zero(W_T); pole(W_T)]);
        fprintf('%2d %-3s = %5.3g: %s\n', link, param, value, ...
            sprintf('%.3g ', sort(roots_abs(roots_abs > 0))));
    end
    for ax = axes_k
        set(ax, 'XScale', 'log');
        grid(ax, 'on');
        title(ax, sprintf('%d. %s', link, names{link}));
        xlabel(ax, '\omega, рад/с');
        legend(ax, 'Location', 'southoutside', 'NumColumns', numel(values), 'FontSize', 8);
    end
    ylabel(axes_k(1), 'L(\omega), дБ');
    ylabel(axes_k(2), '\psi(\omega), град');
end
