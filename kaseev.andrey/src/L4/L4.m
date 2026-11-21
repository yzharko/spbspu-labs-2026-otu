clearvars;
close all;
clc;
set(groot, 'defaultAxesColor', 'w', 'defaultAxesXColor', 'k', ...
    'defaultAxesYColor', 'k', 'defaultTextColor', 'k', ...
    'defaultAxesGridColor', [0.75 0.75 0.75]);

num = 100;
denOpen = [5 0.1 2 2 1];
denClosed = denOpen;
denClosed(end) = denClosed(end) + num;
W = tf(num, denOpen);
H = feedback(W, 1);
pOpen = pole(W);
pClosed = pole(H);

hurwitzOpen = [denOpen(2) denOpen(4) 0 0; ...
    denOpen(1) denOpen(3) denOpen(5) 0; ...
    0 denOpen(2) denOpen(4) 0; ...
    0 denOpen(1) denOpen(3) denOpen(5)];
hurwitzClosed = [denClosed(2) denClosed(4) 0 0; ...
    denClosed(1) denClosed(3) denClosed(5) 0; ...
    0 denClosed(2) denClosed(4) 0; ...
    0 denClosed(1) denClosed(3) denClosed(5)];
dOpen = arrayfun(@(k) det(hurwitzOpen(1:k,1:k)), 1:4);
dClosed = arrayfun(@(k) det(hurwitzClosed(1:k,1:k)), 1:4);

outDir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'figures');
if ~isfolder(outDir)
    mkdir(outDir);
end

f = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1050 720]);
nyquist(W, logspace(-3, 3, 2000));
grid on;
mainTitle = title('Годограф Найквиста разомкнутой системы, вариант 4');
mainTitle.Color = 'k';
set(findall(f, 'Type', 'axes'), 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
set(findall(f, 'Type', 'text'), 'Color', 'k');
exportgraphics(f, fullfile(outDir, 'nyquist.png'), 'Resolution', 180);
close(f);

omegaN = linspace(-10000, 10000, 1000001);
sN = 1i*omegaN;
response = num ./ polyval(denOpen, sN);
argOnePlusW = unwrap(angle(1 + response));
winding = round((argOnePlusW(end)-argOnePlusW(1))/(2*pi));
P = sum(real(pOpen) > 0);
Z = sum(real(pClosed) > 0);
assert(winding == 0 && P == 2 && Z == 2);

omegaM = linspace(0, 6, 3001);
dM = polyval(denClosed, 1i*omegaM);
f = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 620]);
tl = tiledlayout(1, 2, 'Padding', 'compact');
nexttile;
plot(real(dM), imag(dM), 'b-', 'LineWidth', 1.7);
hold on;
plot(101, 0, 'ko', 'MarkerFaceColor', 'k');
grid on;
xlim([90 180]);
ylim([0 4]);
xlabel('Re Dз(jω)');
ylabel('Im Dз(jω)');
title('Начальный участок');
nexttile;
plot(real(dM), imag(dM), 'b-', 'LineWidth', 1.7);
hold on;
plot(101, 0, 'ko', 'MarkerFaceColor', 'k');
grid on;
xlim([0 6500]);
ylim([-12 5]);
xlabel('Re Dз(jω)');
ylabel('Im Dз(jω)');
title('Весь диапазон 0–6 рад/с');
mainTitle = title(tl, 'Кривая Михайлова замкнутой системы');
mainTitle.Color = 'k';
exportgraphics(f, fullfile(outDir, 'mikhailov.png'), 'Resolution', 180);
close(f);

tOpen = 0:0.01:20;
tClosed = 0:0.005:8;
yOpen = step(W, tOpen);
yClosed = step(H, tClosed);
f = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 780]);
tl = tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile;
plot(tOpen, yOpen, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Время, с');
ylabel('yр(t)');
title('Разомкнутая система: расходящиеся колебания');
nexttile;
plot(tClosed, yClosed, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Время, с');
ylabel('yз(t)');
title('Замкнутая система: расходящиеся колебания');
mainTitle = title(tl, 'Реакции на единичное ступенчатое воздействие');
mainTitle.Color = 'k';
exportgraphics(f, fullfile(outDir, 'step_responses.png'), 'Resolution', 180);
close(f);

fid = fopen(fullfile(outDir, 'results.txt'), 'w');
fprintf(fid, 'Variant 4: W(s)=100/(5s^4+0.1s^3+2s^2+2s+1)\n');
fprintf(fid, 'Closed denominator: 5s^4+0.1s^3+2s^2+2s+101\n');
fprintf(fid, 'Open poles:\n');
fprintf(fid, '  %.9f %+.9fi\n', [real(pOpen).'; imag(pOpen).']);
fprintf(fid, 'Closed poles:\n');
fprintf(fid, '  %.9f %+.9fi\n', [real(pClosed).'; imag(pClosed).']);
fprintf(fid, 'Open Hurwitz minors: %.9f %.9f %.9f %.9f\n', dOpen);
fprintf(fid, 'Closed Hurwitz minors: %.9f %.9f %.9f %.9f\n', dClosed);
fprintf(fid, 'Nyquist: P=%d, winding=%d, Z=%d\n', P, winding, Z);
fprintf(fid, 'Mikhailov Im(D)=0 at positive omega=%.9f rad/s\n', sqrt(20));
fclose(fid);

assert(all(dOpen(2:4) < 0));
assert(all(dClosed(2:4) < 0));
disp(fileread(fullfile(outDir, 'results.txt')));
