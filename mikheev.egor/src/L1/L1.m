clear; clc; close all;
K0 = 2.8; T0 = 0.6; T10 = 0.6; T20 = 0.5; T1a = 1.10;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
figs = fullfile(root, 'materials', 'figs');
data_dir = fullfile(root, 'materials', 'results');
if ~exist(figs, 'dir'), mkdir(figs); end
if ~exist(data_dir, 'dir'), mkdir(data_dir); end
syms s t K
s_tf = tf('s');
Ta = (T1a + sqrt(T1a^2 - 4*T20^2))/2;
Tb = (T1a - sqrt(T1a^2 - 4*T20^2))/2;
alpha = T10/(2*T20^2);
omega = sqrt(1/T20^2 - alpha^2);
names = {'1. Интегрирующее', '2. Дифференцирующее', '3. Усилительное', ...
         '4. Апериодическое 1-го порядка', '5. Апериодическое 2-го порядка', ...
         '6. Колебательное', '7. Консервативное', '8. Инерционное интегрирующее', ...
         '9. Реальное дифференцирующее', '10. Форсирующее 1-го порядка', ...
         '11. Изодромное'};
Ws = {K/s, K*s, K, K/(T0*s + 1), K/(T20^2*s^2 + T1a*s + 1), ...
      K/(T20^2*s^2 + T10*s + 1), K/(T0*s^2 + 1), K/(s*(T0*s + 1)), ...
      K*s/(T0*s + 1), K*(T0*s + 1), K*(T0*s + 1)/s};
Wf = {@(k) k/s_tf, @(k) k*s_tf, @(k) tf(k), @(k) k/(T0*s_tf + 1), ...
      @(k) k/(T20^2*s_tf^2 + T1a*s_tf + 1), ...
      @(k) k/(T20^2*s_tf^2 + T10*s_tf + 1), ...
      @(k) k/(T0*s_tf^2 + 1), @(k) k/(s_tf*(T0*s_tf + 1)), ...
      @(k) k*s_tf/(T0*s_tf + 1), @(k) k*(T0*s_tf + 1), ...
      @(k) k*(T0*s_tf + 1)/s_tf};
Tends = [10, 10, 10, 6*T0, 7*max(Ta, Tb), 6*2*pi/omega, ...
         4*2*pi*sqrt(T0), 10, 6*T0, 10, 10];
bad = [2 10];
note_h = {'', 'h(t) = K\delta(t)', '', '', '', '', '', '', '', ...
          'h(t) = K + KT\delta(t)', ''};
note_w = {'', 'w(t) = K\delta''(t)', '+ K\delta(t)', '', '', '', '', ...
          '', '+ (K/T)\delta(t)', 'w(t) = K\delta(t) + KT\delta''(t)', ...
          '+ KT\delta(t)'};
models = repmat(struct('title', '', 'Ws', [], 'Wf', [], 'Tend', 10, ...
                       'hn', '', 'wn', '', 'h', [], 'w', []), 1, 11);
for i = 1:numel(names)
    models(i).title = names{i};
    models(i).Ws = Ws{i};
    models(i).Wf = Wf{i};
    models(i).Tend = Tends(i);
    models(i).hn = note_h{i};
    models(i).wn = note_w{i};
end
fid = fopen(fullfile(data_dir, 'analytic.txt'), 'w');
for i = 1:11
    models(i).h = simplify(ilaplace(models(i).Ws/s, s, t));
    models(i).w = simplify(ilaplace(models(i).Ws, s, t));
    fprintf(fid, '%d|%s|%s|%s\n', i, char(models(i).Ws), ...
            char(models(i).h), char(models(i).w));
end
fclose(fid);
H = cell(1, 11);
W = cell(1, 11);
for i = 1:11
    [tt, y] = from_sym(models(i).h, K0, models(i).Tend);
    H{i} = {tt, y};
    [tt, y] = from_sym(models(i).w, K0, models(i).Tend);
    W{i} = {tt, y};
end
s = struct('fname', 'fig01_h.png', 'ylab', 'h(t)', 'color', [0.05 0.25 0.60], ...
           'notes', {note_h}, 'ylims', {cell(1, 11)}, 'names', {names}, 'dir', figs);
grid_axes(H, s);
s = struct('fname', 'fig02_w.png', 'ylab', 'w(t)', 'color', [0.60 0.10 0.05], ...
           'notes', {note_w}, 'ylims', {cell(1, 11)}, 'names', {names}, 'dir', figs);
grid_axes(W, s);
for mode = 1:2
    C = cell(1, 11);
    L = cell(1, 11);
    YL = cell(1, 11);
    for i = 1:11
        tt = linspace(0, Tends(i), 1500);
        if ismember(i, bad)
            if mode == 1
                if i == 2
                    y = zeros(size(tt));
                    y(1) = K0/(tt(2) - tt(1));
                    L{i} = 'K\delta(t)';
                else
                    y = K0*ones(size(tt));
                    y(1) = K0 + K0*T0/(tt(2) - tt(1));
                    L{i} = 'K + KT\delta(t)';
                end
                YL{i} = [-0.3 2*K0];
            else
                y = zeros(size(tt));
                L{i} = 'проверка неприменима';
            end
        elseif mode == 1
            y = step(Wf{i}(K0), tt);
        else
            y = impulse(Wf{i}(K0), tt);
        end
        C{i} = {tt, y};
    end
    s = struct('fname', 'fig03_step_check.png', 'ylab', 'y', ...
               'color', [0.05 0.25 0.60], 'notes', {L}, 'ylims', {YL}, ...
               'names', {names}, 'dir', figs);
    if mode == 2
        s.fname = 'fig04_impulse_check.png';
        s.color = [0.60 0.10 0.05];
    end
    grid_axes(C, s);
end
Kvals = K0*[0.5 1 1.5 2 3];
s = struct('fname', 'fig05_h_K_link4.png', 'kind', 'step', ...
           'ttl', 'Звено №4: h(t) при изменении K', 'ylab', 'h(t)', ...
           'loc', 'southeast', 'ylims', [], 'xlims', [], ...
           'tt', linspace(0, 6*T0, 1000), 'fmt', 'K = %g', 'dir', figs);
family_axes(Kvals, @(v) tf(v, [T0 1]), s);
s = struct('fname', 'fig06_h_K_link6.png', 'kind', 'step', ...
           'ttl', 'Звено №6: h(t) при изменении K', 'ylab', 'h(t)', ...
           'loc', 'southeast', 'ylims', [], 'xlims', [], ...
           'tt', linspace(0, 6*2*pi/omega, 2000), 'fmt', 'K = %g', 'dir', figs);
family_axes(Kvals, @(v) tf(v, [T20^2 T10 1]), s);
Tvals = T0*[1/5 1/2 1 2 5];
Tmax = 6*max(Tvals);
s = struct('fname', 'fig07_h_T_link4.png', 'kind', 'step', ...
           'ttl', 'Звено №4: h(t) при изменении T', 'ylab', 'h(t)', ...
           'loc', 'southeast', 'ylims', [], 'xlims', [], ...
           'tt', linspace(0, Tmax, 3000), 'fmt', 'T = %g с', 'dir', figs);
family_axes(Tvals, @(v) tf(K0, [v 1]), s);
s = struct('fname', 'fig08_w_T_link4.png', 'kind', 'impulse', ...
           'ttl', 'Звено №4: w(t) при изменении T', 'ylab', 'w(t)', ...
           'loc', 'northeast', 'ylims', [0 2*K0/T0], 'xlims', [], ...
           'tt', linspace(0, Tmax, 3000), 'fmt', 'T = %g с', 'dir', figs);
family_axes(Tvals, @(v) tf(K0, [v 1]), s);
T1vals = T10*[1 0.5 0.25 0.1];
T1max = max(60, 2*6*T20^2/min(T1vals));
xlims = [0 min(T1max, 6*2*pi/omega)];
s = struct('fname', 'fig09_h_T1_link6.png', 'kind', 'step', ...
           'ttl', 'Звено №6: h(t) при уменьшении T1', 'ylab', 'h(t)', ...
           'loc', 'southeast', 'ylims', [], 'xlims', xlims, ...
           'tt', linspace(0, T1max, 4000), 'fmt', 'T1 = %g с', 'dir', figs);
family_axes(T1vals, @(v) tf(K0, [T20^2 v 1]), s);
s = struct('fname', 'fig10_w_T1_link6.png', 'kind', 'impulse', ...
           'ttl', 'Звено №6: w(t) при уменьшении T1', 'ylab', 'w(t)', ...
           'loc', 'northeast', 'ylims', [], 'xlims', xlims, ...
           'tt', linspace(0, T1max, 4000), 'fmt', 'T1 = %g с', 'dir', figs);
family_axes(T1vals, @(v) tf(K0, [T20^2 v 1]), s);
fprintf('\nГотово: рисунки в %s, таблицы в %s\n', figs, data_dir);
function f = canvas(w, h)
    f = figure('Visible', 'off', 'Units', 'centimeters', ...
               'Position', [1 1 w h], 'Color', 'white');
    set(f, 'DefaultAxesFontName', 'Liberation Serif', ...
           'DefaultTextFontName', 'Liberation Serif', ...
           'DefaultAxesFontSize', 8);
end
function [tt, y] = from_sym(expr, k0, tmax)
    ts = sym('t');
    e = subs(expr, sym('K'), k0);
    e = subs(e, {heaviside(ts), dirac(ts), dirac(1, ts)}, {1, 0, 0});
    tt = linspace(0, tmax, 2000).';
    fn = matlabFunction(e, 'Vars', ts);
    y = fn(tt) .* ones(size(tt));
end
function grid_axes(curves, spec)
    f = canvas(17, 15);
    for i = 1:numel(curves)
        subplot(4, 3, i); hold on; grid on; box on;
        tt = curves{i}{1}; y = curves{i}{2};
        plot(tt, y, 'LineWidth', 1.2, 'Color', spec.color);
        if ~isempty(spec.notes{i})
            text(0.03, 0.93, spec.notes{i}, 'Units', 'normalized', ...
                 'FontSize', 6.5, 'Interpreter', 'tex');
        end
        if ~isempty(spec.ylims{i}), ylim(spec.ylims{i}); end
        title(spec.names{i}, 'FontSize', 7);
        xlabel('t, с'); ylabel(spec.ylab); xlim([0 tt(end)]);
        if max(abs(y)) < 1e-9, ylim([-0.5 0.5]); end
        hold off;
    end
    exportgraphics(f, fullfile(spec.dir, spec.fname), 'Resolution', 300);
end
function family_axes(values, make_sys, spec)
    f = canvas(16, 10); hold on; grid on; box on;
    for v = values
        sys = make_sys(v);
        if strcmp(spec.kind, 'step')
            [y, tv] = step(sys, spec.tt);
        else
            [y, tv] = impulse(sys, spec.tt);
        end
        plot(tv, y, 'LineWidth', 1.3, 'DisplayName', sprintf(spec.fmt, v));
    end
    if ~isempty(spec.ylims), ylim(spec.ylims); end
    if ~isempty(spec.xlims), xlim(spec.xlims); end
    title(spec.ttl, 'FontName', 'Liberation Serif');
    xlabel('t, с'); ylabel(spec.ylab);
    legend('Location', spec.loc); legend boxoff;
    exportgraphics(f, fullfile(spec.dir, spec.fname), 'Resolution', 300);
end
