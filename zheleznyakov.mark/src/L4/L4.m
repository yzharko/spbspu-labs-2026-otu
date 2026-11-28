clear; clc; close all;

K = 10;
den = [2 2 1 1];
num = K;
n = numel(den) - 1;
tol = 1e-8;
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
figdir = fullfile(root, 'materials', 'figs');
resdir = fullfile(root, 'materials', 'results');
if ~exist(figdir, 'dir'), mkdir(figdir); end
if ~exist(resdir, 'dir'), mkdir(resdir); end

W = tf(num, den);
Wc = feedback(W, 1);
denc = den; denc(end) = denc(end) + K;
p = sortpoles(pole(W));
pc = sortpoles(pole(Wc));
[H, d] = hurwitz(den);
[Hc, dc] = hurwitz(denc);

aRows = NaN(2, 6);
aRows(1,1:numel(den)) = den;
aRows(2,1:numel(denc)) = denc;
dRows = NaN(2, 5);
dRows(1,1:numel(d)) = d;
dRows(2,1:numel(dc)) = dc;
T = table({'open'; 'closed'}, [n; n], ...
    aRows(:,1), aRows(:,2), aRows(:,3), aRows(:,4), aRows(:,5), aRows(:,6), ...
    dRows(:,1), dRows(:,2), dRows(:,3), dRows(:,4), dRows(:,5), ...
    [double(all(d > tol)); double(all(dc > tol))], ...
    [sum(real(p) > tol); sum(real(pc) > tol)], ...
    'VariableNames', {'system','order','a0','a1','a2','a3','a4','a5', ...
    'd1','d2','d3','d4','d5','stable','rhp'});
writetable(T, fullfile(resdir, 'hurwitz.csv'));

T = table([repmat({'open'}, numel(p), 1); repmat({'closed'}, numel(pc), 1)], ...
    [(1:numel(p))'; (1:numel(pc))'], [real(p); real(pc)], [imag(p); imag(pc)], ...
    'VariableNames', {'system','index','re','im'});
writetable(T, fullfile(resdir, 'roots.csv'));

wNy = [0; logspace(-3, 3, 4000).'];
Hny = freqresp(W, wNy);
Hny = Hny(:);
dist = abs(Hny + 1);
[minDist, kMin] = min(dist);
wStar = wNy(kMin);
P = sum(real(p) > tol);
axisPoles = sum(abs(real(p)) <= tol);
N = winding(num, den);
Z = P + N;
T = table({'open'}, P, N, Z, axisPoles, wStar, minDist, ...
    'VariableNames', {'system','P','N','Z','poles_on_axis','omega_star','min_dist'});
writetable(T, fullfile(resdir, 'nyquist.csv'));

wM = [0; logspace(-3, log10(200), 4000).'];
Dc = polyval(denc, 1j*wM);
argPi = sum(angle(Dc(2:end)./Dc(1:end-1)))/pi;
T = table({'closed'}, n, argPi, sum(real(pc) > tol), double(all(dc > tol)), ...
    double(any(abs(real(pc)) <= tol)), ...
    'VariableNames', {'system','order','arg_pi','rhp','stable','through_origin'});
writetable(T, fullfile(resdir, 'mikhailov.csv'));

tag = {'open'; 'closed'};
tt = cell(2,1); yy = cell(2,1); dcG = zeros(2,1); yEnd = zeros(2,1);
yMax = zeros(2,1); divg = zeros(2,1); oscl = zeros(2,1); tEnd = zeros(2,1);
for k = 1:2
    if k == 1
        pk = p; Wk = W;
    else
        pk = pc; Wk = Wc;
    end
    if max(real(pk)) > tol
        tEnd(k) = 5/max(real(pk));
    else
        tEnd(k) = 3*2*pi/max(abs(imag(pk)));
    end
    tt{k} = linspace(0, tEnd(k), 4000).';
    yy{k} = step(Wk, tt{k});
    dcG(k) = dcgain(Wk);
    yEnd(k) = yy{k}(end);
    yMax(k) = max(abs(yy{k}));
    divg(k) = double(any(real(pk) > tol));
    oscl(k) = double(~divg(k) && any(abs(real(pk)) <= tol));
end
T = table(tag, dcG, yEnd, yMax, divg, oscl, tEnd, ...
    'VariableNames', {'system','dc_gain','y_end','max_abs','diverges','oscillates','t_end'});
writetable(T, fullfile(resdir, 'step.csv'));

figure('Color', 'w');
nyquist(W);
grid off; set(gca, 'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.2);
axis equal; hold on;
xline(0, ':'); yline(0, ':');
plot(-1, 0, 'ro', 'MarkerSize', 8, 'LineWidth', 1.4);
text(-1.2, 0.5, '(-1; 0)');
xlabel('Re W(j\omega)'); ylabel('Im W(j\omega)');
title('Годограф Найквиста разомкнутой системы', 'FontWeight', 'normal');
xlim([-14 14]); ylim([-14 14]);
exportgraphics(gcf, fullfile(figdir, 'fig01_nyquist_full.png'), 'Resolution', 300);
close(gcf);

figure('Color', 'w');
nyquist(W);
grid off; set(gca, 'XGrid', 'on', 'YGrid', 'on', 'GridAlpha', 0.2);
axis equal; hold on;
xline(0, ':'); yline(0, ':');
plot(-1, 0, 'ro', 'MarkerSize', 8, 'LineWidth', 1.4);
plot(real(Hny(kMin)), imag(Hny(kMin)), 'go', 'MarkerSize', 5, ...
    'MarkerFaceColor', 'g');
text(-1.2, 0.35, '(-1; 0)');
xlabel('Re W(j\omega)'); ylabel('Im W(j\omega)');
title('Окрестность точки (-1; 0)', 'FontWeight', 'normal');
xlim([-4.5 1.5]); ylim([-3.5 3.5]);
exportgraphics(gcf, fullfile(figdir, 'fig02_nyquist_zoom.png'), 'Resolution', 300);
close(gcf);

wFull = linspace(0, 200, 20000);
Dfull = polyval(denc, 1j*wFull);
figure('Color', 'w');
plot(real(Dfull), imag(Dfull), 'LineWidth', 1.5);
hold on; grid on; box on;
plot(real(Dfull(1)), imag(Dfull(1)), 'o', 'MarkerFaceColor', 'b');
plot(real(Dfull(end)), imag(Dfull(end)), 's', 'MarkerFaceColor', 'r');
xlabel('Re D_c(j\omega)'); ylabel('Im D_c(j\omega)');
title('Кривая Михайлова замкнутой системы');
exportgraphics(gcf, fullfile(figdir, 'fig03_mikhailov_full.png'), 'Resolution', 300);
close(gcf);

wz = linspace(0, 2.7, 800);
Dz = polyval(denc, 1j*wz);
kRe = find(diff(sign(real(Dz))) ~= 0, 1);
kIm = find(diff(sign(imag(Dz))) ~= 0, 1);
figure('Color', 'w');
plot(real(Dz), imag(Dz), 'LineWidth', 1.5);
hold on; grid on; box on;
plot(real(Dz(1)), imag(Dz(1)), 'o', 'MarkerFaceColor', 'b');
plot(real(Dz(kRe)), imag(Dz(kRe)), 'o', 'MarkerFaceColor', 'b');
plot(real(Dz(kIm)), imag(Dz(kIm)), 's', 'MarkerFaceColor', 'r');
text(real(Dz(kRe)) + 0.3, imag(Dz(kRe)) - 1, ...
    sprintf('\\omega = %.3g', wz(kRe)));
text(real(Dz(kIm)) + 0.8, imag(Dz(kIm)) + 1, ...
    sprintf('\\omega = %.3g', wz(kIm)));
xlabel('Re D_c(j\omega)'); ylabel('Im D_c(j\omega)');
title('Начальный участок кривой Михайлова');
xlim([-3 12.5]); ylim([-26 2]);
exportgraphics(gcf, fullfile(figdir, 'fig04_mikhailov_zoom.png'), 'Resolution', 300);
close(gcf);

figure('Color', 'w');
subplot(1, 2, 1);
step(W, tt{1}); grid on; box on;
title('Разомкнутая система');
xlabel('t, с'); ylabel('y(t)');
subplot(1, 2, 2);
step(Wc, tt{2}); grid on; box on;
title('Замкнутая система');
xlabel('t, с'); ylabel('y(t)');
exportgraphics(gcf, fullfile(figdir, 'fig05_step.png'), 'Resolution', 300);
close(gcf);

tCmp = min(tEnd);
figure('Color', 'w');
step(W, Wc, tCmp); grid on; box on;
legend('разомкнутая', 'замкнутая', 'Location', 'northwest');
xlabel('t, с'); ylabel('y(t)');
title('Реакции на ступенчатое воздействие');
exportgraphics(gcf, fullfile(figdir, 'fig06_step_compare.png'), 'Resolution', 300);
close(gcf);

function p = sortpoles(p)
[~, k] = sortrows([-real(p), -imag(p)]);
p = p(k);
end

function [H, d] = hurwitz(a)
n = numel(a) - 1;
H = zeros(n);
for i = 1:n
    for j = 1:n
        k = 2*j - i + 1;
        if k >= 1 && k <= n + 1
            H(i,j) = a(k);
        end
    end
end
d = zeros(1, n);
for k = 1:n
    d(k) = det(H(1:k, 1:k));
end
d(abs(d) < 1e-12) = 0;
end

function N = winding(num, den)
p = roots(den);
y = sort(imag(p(abs(real(p)) <= 1e-8)));
s = [];
prev = -1e4;
for k = 1:numel(y)
    s = [s; 1j*linspace(prev, y(k) - 1e-3, 800).'];
    s = [s; 1j*y(k) + 1e-3*exp(1j*linspace(-pi/2, pi/2, 120).')];
    prev = y(k) + 1e-3;
end
s = [s; 1j*linspace(prev, 1e4, 800).'];
s = [s; 1e4*exp(1j*linspace(pi/2, -pi/2, 160).')];
denc = den;
denc(end) = denc(end) + num;
f = polyval(denc, s)./polyval(den, s);
N = round(-sum(angle(f(2:end)./f(1:end-1)))/(2*pi));
end
