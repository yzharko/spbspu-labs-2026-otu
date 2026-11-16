clear; clc; close all;
K0 = 1.5; T0 = 0.4; T10 = 0.6; T20 = 0.8; T1a = 1.76;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
figdir = fullfile(root, 'materials', 'figs');
resdir = fullfile(root, 'materials', 'results');
if ~exist(figdir, 'dir'), mkdir(figdir); end
if ~exist(resdir, 'dir'), mkdir(resdir); end
syms s t K
s_tf = tf('s');
Ta = (T1a + sqrt(T1a^2 - 4*T20^2))/2;
Tb = (T1a - sqrt(T1a^2 - 4*T20^2))/2;
alpha = T10/(2*T20^2);
omega = sqrt(1/T20^2 - alpha^2);
titles = {'1. Интегрирующее', '2. Дифференцирующее', '3. Усилительное', ...
          '4. Апериодическое 1-го порядка', '5. Апериодическое 2-го порядка', ...
          '6. Колебательное', '7. Консервативное', '8. Инерционное интегрирующее', ...
          '9. Реальное дифференцирующее', '10. Форсирующее 1-го порядка', ...
          '11. Изодромное'};
titles_nr = numel(titles);
Wsym = {K/s, K*s, K, K/(T0*s + 1), K/(T20^2*s^2 + T1a*s + 1), ...
        K/(T20^2*s^2 + T10*s + 1), K/(T0*s^2 + 1), K/(s*(T0*s + 1)), ...
        K*s/(T0*s + 1), K*(T0*s + 1), K*(T0*s + 1)/s};
Wnum = {@(k) k/s_tf, @(k) k*s_tf, @(k) tf(k), @(k) k/(T0*s_tf + 1), ...
        @(k) k/(T20^2*s_tf^2 + T1a*s_tf + 1), ...
        @(k) k/(T20^2*s_tf^2 + T10*s_tf + 1), ...
        @(k) k/(T0*s_tf^2 + 1), @(k) k/(s_tf*(T0*s_tf + 1)), ...
        @(k) k*s_tf/(T0*s_tf + 1), @(k) k*(T0*s_tf + 1), ...
        @(k) k*(T0*s_tf + 1)/s_tf};
Tend = [10, 10, 10, 6*T0, 7*max(Ta, Tb), 6*2*pi/omega, ...
        4*2*pi*sqrt(T0), 10, 6*T0, 10, 10];
improper = [2 10];
hnote = {'', 'h(t) = K\delta(t)', '', '', '', '', '', '', '', ...
         'h(t) = K + KT\delta(t)', ''};
wnote = {'', 'w(t) = K\delta''(t)', '+ K\delta(t)', '', '', '', '', ...
         '', '+ (K/T)\delta(t)', 'w(t) = K\delta(t) + KT\delta''(t)', ...
         '+ KT\delta(t)'};
h_expr = cell(1, titles_nr);
w_expr = cell(1, titles_nr);
fid = fopen(fullfile(resdir, 'analytic.txt'), 'w');
for i = 1:numel(titles)
    h_expr{i} = simplify(ilaplace(Wsym{i}/s, s, t));
    w_expr{i} = simplify(ilaplace(Wsym{i}, s, t));
    fprintf(fid, '%d|%s|%s|%s\n', i, char(Wsym{i}), ...
            char(h_expr{i}), char(w_expr{i}));
end
fclose(fid);
H = cell(1, titles_nr);
W = cell(1, titles_nr);
for i = 1:titles_nr
    [tt, y] = time_curve(h_expr{i}, K0, Tend(i));
    H{i} = {tt, y};
    [tt, y] = time_curve(w_expr{i}, K0, Tend(i));
    W{i} = {tt, y};
end
s = struct('fname', 'fig01_h.png', 'ylab', 'h(t)', 'color', [0.05 0.25 0.60], ...
           'notes', {hnote}, 'ylims', {cell(1, titles_nr)}, 'titles', {titles}, 'dir', figdir);
draw_tiles(H, s);
s = struct('fname', 'fig02_w.png', 'ylab', 'w(t)', 'color', [0.60 0.10 0.05], ...
           'notes', {wnote}, 'ylims', {cell(1, titles_nr)}, 'titles', {titles}, 'dir', figdir);
draw_tiles(W, s);
for mode = 1:2
    plot_data = cell(1, titles_nr);
    notes = cell(1, titles_nr);
    y_limits = cell(1, titles_nr);
    for i = 1:titles_nr
        tt = linspace(0, Tend(i), 1500);
        if ismember(i, improper)
            if mode == 1
                if i == 2
                    y = zeros(size(tt));
                    y(1) = K0/(tt(2) - tt(1));
                    notes{i} = 'K\delta(t)';
                else
                    y = K0*ones(size(tt));
                    y(1) = K0 + K0*T0/(tt(2) - tt(1));
                    notes{i} = 'K + KT\delta(t)';
                end
                y_limits{i} = [-0.3 2*K0];
            else
                y = zeros(size(tt));
                notes{i} = 'проверка неприменима';
            end
        else
            if mode == 1
                y = step(Wnum{i}(K0), tt);
            else
                y = impulse(Wnum{i}(K0), tt);
            end
        end
        plot_data{i} = {tt, y};
    end
    s = struct('fname', 'fig03_step_check.png', 'ylab', 'y', ...
               'color', [0.05 0.25 0.60], 'notes', {notes}, 'ylims', {y_limits}, ...
               'titles', {titles}, 'dir', figdir);
    if mode == 2
        s.fname = 'fig04_impulse_check.png';
        s.color = [0.60 0.10 0.05];
    end
    draw_tiles(plot_data, s);
end
Kset = K0*[0.5 1 1.5 2 3];
for mode = 1:2
    f = blank_figure(16, 10);
    hold on; grid on; box on;
    for q = 1:numel(Kset)
        if mode == 1
            Wq = tf(Kset(q), [T0 1]);
            tt = linspace(0, 6*T0, 1000);
            ttl = 'Звено №4: h(t) при изменении K';
            fn = 'fig05_h_K_link4.png';
        else
            Wq = tf(Kset(q), [T20^2 T10 1]);
            tt = linspace(0, 6*2*pi/omega, 2000);
            ttl = 'Звено №6: h(t) при изменении K';
            fn = 'fig06_h_K_link6.png';
        end
        [y, tv] = step(Wq, tt);
        plot(tv, y, 'LineWidth', 1.3, 'DisplayName', sprintf('K = %g', Kset(q)));
    end
    title(ttl, 'FontName', 'Times New Roman');
    xlabel('t, с'); ylabel('h(t)');
    legend('Location', 'southeast'); legend boxoff;
    exportgraphics(f, fullfile(figdir, fn), 'Resolution', 300);
end
Tset = T0*[1/5 1/2 1 2 5];
tmaxT = 6*max(Tset);
for mode = 1:2
    f = blank_figure(16, 10);
    hold on; grid on; box on;
    for q = 1:numel(Tset)
        Wq = tf(K0, [Tset(q) 1]);
        tt = linspace(0, tmaxT, 3000);
        if mode == 1
            [y, tv] = step(Wq, tt);
            ttl = 'Звено №4: h(t) при изменении T';
            fn = 'fig07_h_T_link4.png';
            yl = 'h(t)';
            loc = 'southeast';
        else
            [y, tv] = impulse(Wq, tt);
            ttl = 'Звено №4: w(t) при изменении T';
            fn = 'fig08_w_T_link4.png';
            yl = 'w(t)';
            loc = 'northeast';
        end
        plot(tv, y, 'LineWidth', 1.3, 'DisplayName', sprintf('T = %g с', Tset(q)));
    end
    if mode == 2, ylim([0 2*K0/T0]); end
    title(ttl, 'FontName', 'Times New Roman');
    xlabel('t, с'); ylabel(yl);
    legend('Location', loc); legend boxoff;
    exportgraphics(f, fullfile(figdir, fn), 'Resolution', 300);
end
T1set = T10*[1 0.5 0.25 0.1];
tmaxT1 = max(60, 2*6*T20^2/min(T1set));
xlimit = [0 min(tmaxT1, 6*2*pi/omega)];
for mode = 1:2
    f = blank_figure(16, 10);
    hold on; grid on; box on;
    for q = 1:numel(T1set)
        Wq = tf(K0, [T20^2 T1set(q) 1]);
        tt = linspace(0, tmaxT1, 4000);
        if mode == 1
            [y, tv] = step(Wq, tt);
            ttl = 'Звено №6: h(t) при уменьшении T1';
            fn = 'fig09_h_T1_link6.png';
            yl = 'h(t)';
            loc = 'southeast';
        else
            [y, tv] = impulse(Wq, tt);
            ttl = 'Звено №6: w(t) при уменьшении T1';
            fn = 'fig10_w_T1_link6.png';
            yl = 'w(t)';
            loc = 'northeast';
        end
        plot(tv, y, 'LineWidth', 1.3, 'DisplayName', sprintf('T1 = %g с', T1set(q)));
    end
    xlim(xlimit);
    title(ttl, 'FontName', 'Times New Roman');
    xlabel('t, с'); ylabel(yl);
    legend('Location', loc); legend boxoff;
    exportgraphics(f, fullfile(figdir, fn), 'Resolution', 300);
end
fprintf('\nГотово: рисунки в %s, таблицы в %s\n', figdir, resdir);
function f = blank_figure(w, h)
    f = figure('Visible', 'off', 'Units', 'centimeters', ...
               'Position', [1 1 w h], 'Color', 'white');
    set(f, 'DefaultAxesFontName', 'Times New Roman', ...
           'DefaultTextFontName', 'Times New Roman', ...
           'DefaultAxesFontSize', 8);
end
function [tt, y] = time_curve(expr, k0, tmax)
    ts = sym('t');
    e = subs(expr, sym('K'), k0);
    e = subs(e, {heaviside(ts), dirac(ts), dirac(1, ts)}, {1, 0, 0});
    tt = linspace(0, tmax, 2000).';
    f = matlabFunction(e, 'Vars', ts);
    y = f(tt) .* ones(size(tt));
end
function draw_tiles(data, spec)
    f = blank_figure(17, 15);
    tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    for i = 1:numel(data)
        nexttile; hold on; grid on; box on;
        tt = data{i}{1}; y = data{i}{2};
        plot(tt, y, 'LineWidth', 1.2, 'Color', spec.color);
        if ~isempty(spec.notes{i})
            text(0.03, 0.93, spec.notes{i}, 'Units', 'normalized', ...
                 'FontSize', 6.5, 'Interpreter', 'tex');
        end
        if ~isempty(spec.ylims{i}), ylim(spec.ylims{i}); end
        title(spec.titles{i}, 'FontSize', 7);
        xlabel('t, с'); ylabel(spec.ylab); xlim([0 tt(end)]);
        if max(abs(y)) < 1e-9, ylim([-0.5 0.5]); end
        hold off;
    end
    exportgraphics(f, fullfile(spec.dir, spec.fname), 'Resolution', 300);
end
