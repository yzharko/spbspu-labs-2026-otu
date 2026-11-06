clear; clc; close all;

K = 1.1;
T = 0.2;
T1 = 0.3;
T2 = 0.5;
%у варианта T1 < 2*T2, ближайшее допустимое значение для звена 5: T1 = 2*T2
T1ap = 2 * T2;

rootDir = fileparts(mfilename('fullpath'));
if isempty(rootDir)
    rootDir = pwd;
end
figDir = fullfile(rootDir, 'fig');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

if usejava('desktop')
    set(0, 'DefaultFigureVisible', 'on');
else
    set(0, 'DefaultFigureVisible', 'off');
end
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultTextInterpreter', 'none');
set(0, 'DefaultLegendInterpreter', 'none');
set(0, 'DefaultAxesTickLabelInterpreter', 'none');
set(0, 'DefaultFigureColor', 'w');
set(0, 'DefaultAxesColor', 'w');
set(0, 'DefaultAxesXColor', 'k');
set(0, 'DefaultAxesYColor', 'k');
set(0, 'DefaultTextColor', 'k');
set(0, 'DefaultAxesGridColor', [0.35 0.35 0.35]);
set(0, 'DefaultLegendColor', 'w');
set(0, 'DefaultLegendTextColor', 'k');
set(0, 'DefaultLegendEdgeColor', [0.2 0.2 0.2]);

syms s t
syms ks ts t1s t2s positive
old = [ks ts t1s t2s];

%вход один и тот же для всех звеньев
U_step = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);

fid = fopen(fullfile(rootDir, 'results.txt'), 'w', 'n', 'UTF-8');
logf(fid, 'U_step = %s\n', char(U_step));
logf(fid, 'U_impulse = %s\n', char(U_impulse));
logf(fid, 'K = %g\nT = %g\nT1 = %g\nT2 = %g\nT1ap = %g\n', K, T, T1, T2, T1ap);
logf(fid, 'discriminant = %g\n', T1^2 - 4*T2^2);

links = defineLinks(ks, ts, t1s, t2s, s);
Klist = [0.55 1.1 2.2];

try
for i = 1:numel(links)
    L = links(i);
    p0 = nominalOf(L.mode, K, T, T1, T2, T1ap);
    W = L.W;
    if L.id == 5
        W = subs(L.W, t1s, 2*t2s);
    end
    Y_step = simplify(W * U_step);
    Y_impulse = simplify(W * U_impulse);
    h = simplify(ilaplace(Y_step, s, t));
    w = simplify(ilaplace(Y_impulse, s, t));
    [hReg, h0, h1] = splitImpulses(h, t);
    [wReg, w0, w1] = splitImpulses(w, t);

    logf(fid, '\n=== %d. %s ===\n', L.id, L.name);
    logf(fid, 'W = %s\n', char(simplify(W)));
    logf(fid, 'Y_step = %s\n', char(Y_step));
    logf(fid, 'Y_impulse = %s\n', char(Y_impulse));
    logf(fid, 'h = %s\n', char(h));
    logf(fid, 'w = %s\n', char(w));
    logf(fid, 'h_reg = %s\n', char(hReg));
    logf(fid, 'h_c0 = %s\n', char(h0));
    logf(fid, 'h_c1 = %s\n', char(h1));
    logf(fid, 'w_reg = %s\n', char(wReg));
    logf(fid, 'w_c0 = %s\n', char(w0));
    logf(fid, 'w_c1 = %s\n', char(w1));
    logf(fid, 'h_t>0 = %s\n', char(simplify(subs(subs(h, old, p0), heaviside(t), 1))));
    logf(fid, 'w_t>0 = %s\n', char(simplify(subs(subs(w, old, p0), heaviside(t), 1))));

    tNum = timeGrid(L.tEnd);
    hA = evalSub(hReg, t, tNum, old, p0);
    wA = evalSub(wReg, t, tNum, old, p0);
    checkLink(fid, L.id, tNum, hA, wA, p0);

    [hC, hOk, hMsg] = cstResponse('step', L.W, old, p0, tNum);
    [wC, wOk, wMsg] = cstResponse('impulse', L.W, old, p0, tNum);
    logf(fid, 'check_h ok=%d err=%g msg=%s\n', hOk, curveErr(tNum, hA, hC, hOk), oneLine(hMsg));
    logf(fid, 'check_w ok=%d err=%g msg=%s\n', wOk, curveErr(tNum, wA, wC, wOk), oneLine(wMsg));

    hNote = impulseNote(h0, h1, old, p0);
    wNote = impulseNote(w0, w1, old, p0);
    plotPair(fullfile(figDir, sprintf('hw_%02d.png', L.id)), L.name, tNum, ...
        hA, hC, hOk, wA, wC, wOk, hNote, wNote);

    [idxT, valsT, labT] = varySpec(L.mode, T);
    if L.id == 5
        plotAperiodic2(fullfile(figDir, 'var_05.png'), L.name, Klist, T2, L.tEndVar);
    else
        plotVariations(fullfile(figDir, sprintf('var_%02d.png', L.id)), L.name, ...
            hReg, h0, h1, t, old, p0, Klist, idxT, valsT, labT, L.tEndVar);
    end

    if L.id == 4
        logAperiodic(fid, hReg, t, old, p0, T);
    elseif L.id == 6
        logOscillatory(fid, hReg, t, old, p0, T2);
    end
end

logf(fid, '\nDONE\n');
fclose(fid);
catch ME
    fclose(fid);
    rethrow(ME);
end

function L = defineLinks(ks, ts, t1s, t2s, s)
    spec = {
        1, 'Интегрирующее', ks/s, 10, 10, 'K';
        2, 'Дифференцирующее', ks*s, 2, 2, 'K';
        3, 'Усилительное', ks, 2, 2, 'K';
        4, 'Апериодическое 1-го порядка', ks/(ts*s + 1), 3, 8, 'T';
        5, 'Апериодическое 2-го порядка', ks/(t2s^2*s^2 + t1s*s + 1), 8, 20, 'T1ap';
        6, 'Колебательное', ks/(t2s^2*s^2 + t1s*s + 1), 15, 45, 'T1osc';
        7, 'Консервативное', ks/(ts*s^2 + 1), 10, 20, 'T';
        8, 'Инерционное интегрирующее', ks/(s*(ts*s + 1)), 10, 12, 'T';
        9, 'Реальное дифференцирующее', ks*s/(ts*s + 1), 2, 8, 'T';
        10, 'Форсирующее 1-го порядка', ks*(ts*s + 1), 2, 2, 'T';
        11, 'Изодромное', ks*(ts*s + 1)/s, 10, 10, 'T'
    };
    L = struct('id', {}, 'name', {}, 'W', {}, 'tEnd', {}, 'tEndVar', {}, 'mode', {});
    for i = 1:size(spec, 1)
        L(i).id = spec{i, 1};
        L(i).name = spec{i, 2};
        L(i).W = spec{i, 3};
        L(i).tEnd = spec{i, 4};
        L(i).tEndVar = spec{i, 5};
        L(i).mode = spec{i, 6};
    end
end

function p = nominalOf(mode, K, T, T1, T2, T1ap)
    t1 = T1;
    if strcmp(mode, 'T1ap')
        t1 = T1ap;
    end
    p = [K T t1 T2];
end

function [idx, vals, label] = varySpec(mode, T)
    switch mode
        case 'K'
            idx = [];
            vals = [];
            label = '';
        case 'T'
            idx = 2;
            vals = round(T * [0.2 0.5 1 2 5], 10);
            label = 'T';
        case 'T1osc'
            idx = 3;
            vals = [0.9 0.6 0.3 0.15 0.05];
            label = 'T1';
        case 'T1ap'
            idx = 3;
            vals = [1 1.5 2 4];
            label = 'T1';
        otherwise
            error('неизвестный режим %s', mode);
    end
end

function [reg, c0, c1] = splitImpulses(expr, t)
    d0 = sym('d0');
    d1 = sym('d1');
    e = subs(expr, dirac(0, t), dirac(t));
    e = subs(e, [dirac(t), dirac(1, t)], [d0, d1]);
    c0 = simplify(subs(diff(e, d0), [d0 d1], [0 0]));
    c1 = simplify(subs(diff(e, d1), [d0 d1], [0 0]));
    reg = simplify(subs(e, [d0 d1], [0 0]));
end

function y = evalSub(expr, tSym, tNum, old, new)
    e = subs(expr, [old heaviside(tSym)], [new 1]);
    if ~isreal(e)
        e = simplify(rewrite(e, 'sincos'), 'IgnoreAnalyticConstraints', true);
    end
    y = evalExpr(e, tSym, tNum);
end

function y = evalExpr(expr, tSym, tNum)
    if isempty(symvar(expr))
        y = double(expr) * ones(size(tNum));
    else
        f = matlabFunction(expr, 'Vars', tSym);
        y = f(tNum(:));
        y = reshape(double(y), size(tNum));
    end
    if ~isreal(y)
        if max(abs(imag(y(:)))) < 1e-8
            y = real(y);
        else
            error('комплексные значения на сетке');
        end
    end
end

function tNum = timeGrid(tEnd)
    dt = min(0.002, tEnd / 2500);
    tNum = 0:dt:tEnd;
end

function checkLink(fid, id, tNum, hA, wA, p)
    Kp = p(1);
    Tp = p(2);
    switch id
        case 1
            must(fid, id, 'h', hA, Kp * tNum);
            must(fid, id, 'w', wA, Kp * ones(size(tNum)));
        case 2
            must(fid, id, 'h', hA, zeros(size(tNum)));
            must(fid, id, 'w', wA, zeros(size(tNum)));
        case 3
            must(fid, id, 'h', hA, Kp * ones(size(tNum)));
            must(fid, id, 'w', wA, zeros(size(tNum)));
        case 4
            must(fid, id, 'h', hA, Kp * (1 - exp(-tNum / Tp)));
            must(fid, id, 'w', wA, (Kp / Tp) * exp(-tNum / Tp));
            dh = gradient(hA, tNum);
            mid = 30:numel(tNum)-30;
            logf(fid, 'dw err = %g\n', max(abs(dh(mid) - wA(mid))));
        case 7
            wn = 1 / sqrt(Tp);
            must(fid, id, 'h', hA, Kp * (1 - cos(wn * tNum)));
            must(fid, id, 'w', wA, Kp * wn * sin(wn * tNum));
        case 8
            must(fid, id, 'h', hA, Kp * (tNum - Tp * (1 - exp(-tNum / Tp))));
            must(fid, id, 'w', wA, Kp * (1 - exp(-tNum / Tp)));
        case 9
            must(fid, id, 'h', hA, (Kp / Tp) * exp(-tNum / Tp));
        case 10
            must(fid, id, 'h', hA, Kp * ones(size(tNum)));
        case 11
            must(fid, id, 'h', hA, Kp * (tNum + Tp));
            must(fid, id, 'w', wA, Kp * ones(size(tNum)));
    end
end

function must(fid, id, name, y, yRef)
    e = max(abs(y(:) - yRef(:)));
    logf(fid, 'assert %d %s err=%g\n', id, name, e);
    if e > 1e-5
        error('звено %d, %s, погрешность %g', id, name, e);
    end
end

function [y, ok, msg] = cstResponse(kind, W, old, new, tNum)
    y = [];
    ok = false;
    msg = '';
    try
        Wtf = symToTf(W, old, new);
        if strcmp(kind, 'step')
            y = step(Wtf, tNum(:));
        else
            y = impulse(Wtf, tNum(:));
        end
        y = reshape(double(y), size(tNum));
        if any(~isfinite(y))
            ok = false;
            msg = 'бесконечные значения';
            y = [];
            return;
        end
        ok = true;
    catch ME
        msg = ME.message;
        y = [];
        ok = false;
    end
end

function Wtf = symToTf(W, old, new)
    Wn = simplify(subs(W, old, new));
    [n, d] = numden(Wn);
    Wtf = tf(sym2poly(expand(n)), sym2poly(expand(d)));
end

function e = curveErr(tNum, a, b, ok)
    if ~ok
        e = nan;
        return;
    end
    tNum = tNum(:);
    a = a(:);
    b = b(:);
    mask = tNum >= 0.02;
    e = max(abs(a(mask) - b(mask)));
end

function note = impulseNote(c0, c1, old, p)
    a = double(subs(c0, old, p));
    b = double(subs(c1, old, p));
    parts = {};
    if abs(a) > 1e-10
        parts{end + 1} = sprintf('есть дельта(t) весом %.4g', a); %#ok<AGROW>
    end
    if abs(b) > 1e-10
        parts{end + 1} = sprintf('есть производная дельта(t) весом %.4g', b); %#ok<AGROW>
    end
    if isempty(parts)
        note = '';
    else
        note = strjoin(parts, ', ');
    end
end

function plotPair(path, name, tNum, hA, hC, hOk, wA, wC, wOk, hNote, wNote)
    fig = figure('Color', 'w', 'Position', [80 80 1000 430]);
    subplot(1, 2, 1);
    drawCurves(tNum, hA, hC, hOk, 'step', hNote);
    grid on; box on;
    xlabel('t, с');
    ylabel('h(t)');
    title('переходная характеристика');
    subplot(1, 2, 2);
    drawCurves(tNum, wA, wC, wOk, 'impulse', wNote);
    grid on; box on;
    xlabel('t, с');
    ylabel('w(t)');
    title('импульсная характеристика');
    sgtitle(name, 'Interpreter', 'none', 'FontName', 'Times New Roman', 'FontSize', 14);
    styleFig(fig);
    exportgraphics(fig, path, 'Resolution', 160);
    if ~usejava('desktop')
        close(fig);
    end
end

function drawCurves(tNum, yReg, yCst, ok, cstName, note)
    hasReg = max(abs(yReg)) > 1e-8;
    if ~hasReg
        xlim([tNum(1) tNum(end)]);
        ylim([0 1]);
        if isempty(note)
            note = 'регулярная часть равна 0';
        end
        text(0.5, 0.62, note, 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
        if ok
            extra = 'CST даёт только регулярную часть, без импульса';
        else
            extra = 'step/impulse эту реакцию не строит';
        end
        text(0.5, 0.40, extra, 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
        return;
    end
    hold on;
    plot(tNum, yReg, 'b', 'LineWidth', 1.6, 'DisplayName', 'аналитически');
    if ok
        plot(tNum, yCst, 'r--', 'LineWidth', 1.3, 'DisplayName', cstName);
    end
    if ~isempty(note)
        text(0.02, 0.95, note, 'Units', 'normalized', ...
            'VerticalAlignment', 'top', 'FontSize', 9);
    end
    legend('Location', 'best');
    hold off;
end

function plotVariations(path, name, hReg, h0, h1, tSym, old, base, Klist, idxT, valsT, labT, tEnd)
    fig = figure('Color', 'w', 'Position', [70 70 1000 430]);
    if isempty(idxT)
        drawFamily(hReg, h0, h1, tSym, old, base, 1, Klist, 'K', tEnd, 'h(t)');
        title([name ', влияние K'], 'Interpreter', 'none');
        styleFig(fig);
    else
        subplot(1, 2, 1);
        drawFamily(hReg, h0, h1, tSym, old, base, 1, Klist, 'K', tEnd, 'h(t)');
        title('влияние K');
        subplot(1, 2, 2);
        drawFamily(hReg, h0, h1, tSym, old, base, idxT, valsT, labT, tEnd, 'h(t)');
        title(['влияние ' labT]);
        sgtitle(name, 'Interpreter', 'none', 'FontName', 'Times New Roman', 'FontSize', 14);
        styleFig(fig);
    end
    exportgraphics(fig, path, 'Resolution', 160);
    if ~usejava('desktop')
        close(fig);
    end
end

function drawFamily(exprReg, c0, c1, tSym, old, base, idx, vals, label, tEnd, yLabel)
    tNum = timeGrid(tEnd);
    n = numel(vals);
    regs = zeros(n, numel(tNum));
    c0s = zeros(1, n);
    c1s = zeros(1, n);
    for i = 1:n
        p = base;
        p(idx) = vals(i);
        regs(i, :) = evalSub(exprReg, tSym, tNum, old, p);
        c0s(i) = double(subs(c0, old, p));
        c1s(i) = double(subs(c1, old, p));
    end
    spread = max(regs(:)) - min(regs(:));
    hasImpulse = max(abs(c0s)) > 1e-10 || max(abs(c1s)) > 1e-10;
    if max(abs(regs(:))) < 1e-8 || (hasImpulse && spread < 1e-6)
        cats = categorical(compose('%g', vals), compose('%g', vals), 'Ordinal', true);
        if max(abs(c1s)) > 1e-10
            bar(cats, [c0s(:) c1s(:)]);
            legend('вес δ(t)', 'вес δ''(t)', 'Location', 'best');
        else
            bar(cats, c0s);
            ylabel('вес δ(t)');
        end
        xlabel(label);
        grid on; box on;
        return;
    end
    hold on; grid on; box on;
    for i = 1:n
        plot(tNum, regs(i, :), 'LineWidth', 1.35, ...
            'DisplayName', sprintf('%s = %g', label, vals(i)));
    end
    xlabel('t, с');
    ylabel(yLabel);
    legend('Location', 'best');
    hold off;
    if max(abs(c0s)) > 1e-10 || max(abs(c1s)) > 1e-10
        text(0.02, 0.95, 'кроме кривой, в нуле есть импульс', ...
            'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 9);
    end
end

function logAperiodic(fid, hReg, tSym, old, base, T)
    vals = round(T * [0.2 0.5 1 2 5], 10);
    logf(fid, 'METRIC aperiodic header\n');
    t95ref = nan;
    for i = 1:numel(vals)
        p = base;
        p(2) = vals(i);
        tNum = timeGrid(max(8, 8 * vals(i)));
        y = evalSub(hReg, tSym, tNum, old, p);
        t95 = firstCross(tNum, y, 0.95 * p(1));
        if abs(vals(i) - T) < 1e-12
            t95ref = t95;
        end
        logf(fid, 'METRIC ap T=%g t95=%g t95_over_T=%g\n', vals(i), t95, t95 / vals(i));
    end
    for k = [0.55 1.1 2.2]
        p = base;
        p(1) = k;
        tNum = timeGrid(3);
        y = evalSub(hReg, tSym, tNum, old, p);
        t95 = firstCross(tNum, y, 0.95 * k);
        logf(fid, 'METRIC apK K=%g t95=%g\n', k, t95);
    end
    logf(fid, 'METRIC ap t95ref=%g\n', t95ref);
end

function logOscillatory(fid, hReg, tSym, old, base, T2)
    vals = [0.9 0.6 0.3 0.15 0.05];
    Kp = base(1);
    for i = 1:numel(vals)
        p = base;
        p(3) = vals(i);
        zeta = vals(i) / (2 * T2);
        wn = 1 / T2;
        sigA = exp(-pi * zeta / sqrt(1 - zeta^2)) * 100;
        tau = 1 / (zeta * wn);
        tEnd = max(45, 6 * tau);
        tNum = timeGrid(tEnd);
        y = evalSub(hReg, tSym, tNum, old, p);
        sig = (max(y) - Kp) / Kp * 100;
        tset = lastOut(tNum, y, Kp, 0.05 * Kp);
        logf(fid, 'METRIC osc T1=%g zeta=%g sig=%g sigA=%g tset=%g tau=%g wn=%g\n', ...
            vals(i), zeta, sig, sigA, tset, tau, wn);
    end
end

function tCross = firstCross(tNum, y, level)
    idx = find(y >= level, 1, 'first');
    if isempty(idx)
        tCross = inf;
    else
        tCross = tNum(idx);
    end
end

function tSet = lastOut(tNum, y, level, band)
    idx = find(abs(y - level) > band, 1, 'last');
    if isempty(idx)
        tSet = 0;
    elseif idx == numel(tNum)
        tSet = inf;
    else
        tSet = tNum(idx);
    end
end

function plotAperiodic2(path, name, Klist, T2, tEnd)
    t1nom = 2 * T2;
    tNum = timeGrid(tEnd);
    fig = figure('Color', 'w', 'Position', [70 70 1000 430]);
    subplot(1, 2, 1);
    hold on; grid on; box on;
    for k = Klist
        plot(tNum, hSecond(tNum, k, t1nom, T2), 'LineWidth', 1.35, ...
            'DisplayName', sprintf('K = %g', k));
    end
    xlabel('t, с');
    ylabel('h(t)');
    title('влияние K');
    legend('Location', 'best');
    hold off;
    subplot(1, 2, 2);
    hold on; grid on; box on;
    for t1 = [1 1.5 2 4]
        plot(tNum, hSecond(tNum, 1.1, t1, T2), 'LineWidth', 1.35, ...
            'DisplayName', sprintf('T1 = %g', t1));
    end
    xlabel('t, с');
    ylabel('h(t)');
    title('влияние T1');
    legend('Location', 'best');
    hold off;
    sgtitle(name, 'Interpreter', 'none', 'FontName', 'Times New Roman', 'FontSize', 14);
    styleFig(fig);
    exportgraphics(fig, path, 'Resolution', 160);
    if ~usejava('desktop')
        close(fig);
    end
end

function y = hSecond(tNum, K, T1, T2)
    disc = T1^2 - 4 * T2^2;
    if abs(disc) < 1e-9
        a = 1 / T2;
        y = K * (1 - exp(-a * tNum) .* (1 + a * tNum));
    else
        alpha = T1 / (2 * T2^2);
        beta = sqrt(alpha^2 - 1 / T2^2);
        y = K * (1 - exp(-alpha * tNum) .* (cosh(beta * tNum) + (alpha / beta) * sinh(beta * tNum)));
    end
end

function styleFig(fig)
    set(fig, 'Color', 'w');
    axs = findall(fig, 'Type', 'axes');
    for k = 1:numel(axs)
        set(axs(k), 'Color', 'w', 'XColor', 'k', 'YColor', 'k', ...
            'GridColor', [0.45 0.45 0.45], 'GridAlpha', 0.8, ...
            'FontName', 'Times New Roman', 'FontSize', 11);
    end
    txt = findall(fig, 'Type', 'text');
    for k = 1:numel(txt)
        set(txt(k), 'Color', 'k', 'FontName', 'Times New Roman');
    end
    legs = findall(fig, 'Type', 'legend');
    for k = 1:numel(legs)
        set(legs(k), 'Color', 'w', 'TextColor', 'k', ...
            'EdgeColor', [0.3 0.3 0.3], 'FontName', 'Times New Roman');
    end
end

function s = oneLine(msg)
    s = strrep(msg, newline, ' ');
    s = strrep(s, sprintf('\r'), ' ');
end

function logf(fid, fmt, varargin)
    fprintf(fid, fmt, varargin{:});
    fprintf(1, fmt, varargin{:});
end
