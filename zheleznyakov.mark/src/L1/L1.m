clear; clc; close all;

K0 = 1.8; T0 = 0.6; T10 = 0.9; T20 = 1.2; T1a = 2.64;
here = fileparts(mfilename('fullpath')); base = fileparts(fileparts(here));
figdir = fullfile(base, 'materials', 'figs');
resdir = fullfile(base, 'materials', 'results');
if ~exist(figdir, 'dir'), mkdir(figdir); end
if ~exist(resdir, 'dir'), mkdir(resdir); end
fprintf('Вариант: K = %g, T = %g с, T1 = %g с, T2 = %g с, T1a = %g с\n', ...
        K0, T0, T10, T20, T1a);

syms s t K
s_tf = tf('s');
Ta = (T1a + sqrt(T1a^2 - 4*T20^2))/2;
Tb = (T1a - sqrt(T1a^2 - 4*T20^2))/2;
alpha = T10/(2*T20^2);
omega = sqrt(1/T20^2 - alpha^2);

links = {
    struct('name','1. Интегрирующее','W',K/s,'Wc',@(k)k/s_tf, ...
           'improper',false,'tmax',10,'dh','','dw','')
    struct('name','2. Дифференцирующее','W',K*s,'Wc',@(k)k*s_tf, ...
           'improper',true,'tmax',10,'dh','h(t) = K\delta(t)','dw','w(t) = K\delta''(t)')
    struct('name','3. Усилительное','W',K,'Wc',@(k)tf(k), ...
           'improper',false,'tmax',10,'dh','','dw','+ K\delta(t)')
    struct('name','4. Апериодическое 1-го порядка','W',K/(T0*s + 1), ...
           'Wc',@(k)k/(T0*s_tf + 1),'improper',false,'tmax',6*T0,'dh','','dw','')
    struct('name','5. Апериодическое 2-го порядка','W',K/(T20^2*s^2 + T1a*s + 1), ...
           'Wc',@(k)k/(T20^2*s_tf^2 + T1a*s_tf + 1),'improper',false,'tmax',7*max(Ta, Tb),'dh','','dw','')
    struct('name','6. Колебательное','W',K/(T20^2*s^2 + T10*s + 1), ...
           'Wc',@(k)k/(T20^2*s_tf^2 + T10*s_tf + 1),'improper',false,'tmax',6*2*pi/omega,'dh','','dw','')
    struct('name','7. Консервативное','W',K/(T0*s^2 + 1), ...
           'Wc',@(k)k/(T0*s_tf^2 + 1),'improper',false,'tmax',4*2*pi*sqrt(T0),'dh','','dw','')
    struct('name','8. Инерционное интегрирующее','W',K/(s*(T0*s + 1)), ...
           'Wc',@(k)k/(s_tf*(T0*s_tf + 1)),'improper',false,'tmax',10,'dh','','dw','')
    struct('name','9. Реальное дифференцирующее','W',K*s/(T0*s + 1), ...
           'Wc',@(k)k*s_tf/(T0*s_tf + 1),'improper',false,'tmax',6*T0,'dh','','dw','+ (K/T)\delta(t)')
    struct('name','10. Форсирующее 1-го порядка','W',K*(T0*s + 1), ...
           'Wc',@(k)k*(T0*s_tf + 1),'improper',true,'tmax',10, ...
           'dh','h(t) = K + KT\delta(t)','dw','w(t) = K\delta(t) + KT\delta''(t)')
    struct('name','11. Изодромное','W',K*(T0*s + 1)/s, ...
           'Wc',@(k)k*(T0*s_tf + 1)/s_tf,'improper',false,'tmax',10,'dh','','dw','+ KT\delta(t)')
};

N = numel(links);

U_step = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);
fid = fopen(fullfile(resdir, 'analytic.txt'), 'w');
for i = 1:N
    links{i}.h = simplify(ilaplace(links{i}.W*U_step, s, t));
    links{i}.w = simplify(ilaplace(links{i}.W*U_impulse, s, t));
    fprintf('%d. %s\n   W(s) = %s\n   h(t) = %s\n   w(t) = %s\n', ...
            i, links{i}.name, char(links{i}.W), char(links{i}.h), char(links{i}.w));
    fprintf(fid, '%d|%s|%s|%s\n', i, char(links{i}.W), ...
            char(links{i}.h), char(links{i}.w));
end
fclose(fid);

for pass = 1:2
    f = newfig(17, 15);
    tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    for i = 1:N
        nexttile; hold on; grid on; box on;
        if pass == 1
            e = links{i}.h; note = links{i}.dh; col = [0.05 0.25 0.60];
            yl = 'h(t)'; fn = 'fig01_h.png';
        else
            e = links{i}.w; note = links{i}.dw; col = [0.60 0.10 0.05];
            yl = 'w(t)'; fn = 'fig02_w.png';
        end
        [tt, y] = sample_expr(e, K0, links{i}.tmax);
        plot(tt, y, 'LineWidth', 1.2, 'Color', col);
        if ~isempty(note)
            text(0.03, 0.93, note, 'Units', 'normalized', ...
                 'FontSize', 6.5, 'Interpreter', 'tex');
        end
        title(links{i}.name, 'FontSize', 7);
        xlabel('t, с'); ylabel(yl); xlim([0 links{i}.tmax]);
        if max(abs(y)) < 1e-9, ylim([-0.5 0.5]); end
        hold off;
    end
    exportgraphics(f, fullfile(figdir, fn), 'Resolution', 300);
end

for pass = 1:2
    f = newfig(17, 15);
    tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    for i = 1:N
        nexttile; hold on; grid on; box on;
        tt = linspace(0, links{i}.tmax, 1500);
        if pass == 1
            if links{i}.improper
                if i == 2
                    y = zeros(size(tt)); y(1) = K0/(tt(2) - tt(1));
                    note = 'K\delta(t)';
                else
                    y = K0*ones(size(tt)); y(1) = K0 + K0*T0/(tt(2) - tt(1));
                    note = 'K + KT\delta(t)';
                end
                plot(tt, y, 'LineWidth', 1.2);
                text(0.03, 0.93, note, 'Units', 'normalized', ...
                     'FontSize', 6.5, 'Interpreter', 'tex');
                ylim([-0.3 2*K0]);
            else
                y = step(links{i}.Wc(K0), tt);
                plot(tt, y, 'LineWidth', 1.2, 'Color', [0.05 0.25 0.60]);
            end
            fn = 'fig03_step_check.png';
        else
            if links{i}.improper
                plot(tt, zeros(size(tt)), 'LineWidth', 1.0);
                text(0.03, 0.93, 'проверка неприменима', 'Units', 'normalized', ...
                     'FontSize', 6.5, 'Interpreter', 'none');
            else
                y = impulse(links{i}.Wc(K0), tt);
                plot(tt, y, 'LineWidth', 1.2, 'Color', [0.60 0.10 0.05]);
            end
            fn = 'fig04_impulse_check.png';
        end
        title(links{i}.name, 'FontSize', 7);
        xlabel('t, с'); ylabel('y'); xlim([0 links{i}.tmax]);
        hold off;
    end
    exportgraphics(f, fullfile(figdir, fn), 'Resolution', 300);
end

Kfac = [0.5 1 1.5 2 3];
plot_family(figdir, 'fig05_h_K_link4.png', 'Звено №4: h(t) при изменении K', 'h(t)', ...
    linspace(0, 6*T0, 1000), Kfac, @(kf) tf(K0*kf, [T0 1]), 'step', 'K = %g', 'southeast', [], []);
plot_family(figdir, 'fig06_h_K_link6.png', 'Звено №6: h(t) при изменении K', 'h(t)', ...
    linspace(0, 6*2*pi/omega, 2000), Kfac, @(kf) tf(K0*kf, [T20^2 T10 1]), 'step', 'K = %g', 'southeast', [], []);

Tser = T0*[1/5 1/2 1 2 5];
tmaxT = 6*max(Tser);
plot_family(figdir, 'fig07_h_T_link4.png', 'Звено №4: h(t) при изменении T', 'h(t)', ...
    linspace(0, tmaxT, 3000), Tser, @(T) tf(K0, [T 1]), 'step', 'T = %g с', 'southeast', [], []);
plot_family(figdir, 'fig08_w_T_link4.png', 'Звено №4: w(t) при изменении T', 'w(t)', ...
    linspace(0, tmaxT, 3000), Tser, @(T) tf(K0, [T 1]), 'impulse', 'T = %g с', 'northeast', [0 2*K0/T0], []);

T1ser = T10*[1 0.5 0.25 0.1];
tmax6 = max(60, 2*6*T20^2/min(T1ser));
xl = [0 min(tmax6, 6*2*pi/omega)];
plot_family(figdir, 'fig09_h_T1_link6.png', 'Звено №6: h(t) при уменьшении T1', 'h(t)', ...
    linspace(0, tmax6, 4000), T1ser, @(T1) tf(K0, [T20^2 T1 1]), 'step', 'T1 = %g с', 'southeast', [], xl);
plot_family(figdir, 'fig10_w_T1_link6.png', 'Звено №6: w(t) при уменьшении T1', 'w(t)', ...
    linspace(0, tmax6, 4000), T1ser, @(T1) tf(K0, [T20^2 T1 1]), 'impulse', 'T1 = %g с', 'northeast', [], xl);

fprintf('\nГотово: рисунки в %s, таблицы в %s\n', figdir, resdir);

function f = newfig(w, h)
    f = figure('Visible', 'off', 'Units', 'centimeters', ...
               'Position', [1 1 w h], 'Color', 'white');
    set(f, 'DefaultAxesFontName', 'Liberation Serif', ...
           'DefaultTextFontName', 'Liberation Serif', ...
           'DefaultAxesFontSize', 8);
end

function plot_family(figdir, fn, ttl, yl, tt, vals, Wfun, kind, lab, ll, ylims, xlims)
    f = newfig(16, 10); hold on; grid on; box on;
    for v = vals
        W = Wfun(v);
        if strcmp(kind, 'step')
            [y, t] = step(W, tt);
        else
            [y, t] = impulse(W, tt);
        end
        plot(t, y, 'LineWidth', 1.3, 'DisplayName', sprintf(lab, v));
    end
    title(ttl, 'FontName', 'Liberation Serif');
    xlabel('t, с'); ylabel(yl);
    legend('Location', ll); legend boxoff;
    if ~isempty(ylims), ylim(ylims); end
    if ~isempty(xlims), xlim(xlims); end
    exportgraphics(f, fullfile(figdir, fn), 'Resolution', 300);
end

function [tt, y] = sample_expr(expr, k0, tmax)
    ts = sym('t');
    e = subs(expr, sym('K'), k0);
    e = subs(e, {heaviside(ts), dirac(ts), dirac(1, ts)}, {1, 0, 0});
    tt = linspace(0, tmax, 2000).';
    f = matlabFunction(e, 'Vars', ts);
    y = f(tt) .* ones(size(tt));
end
