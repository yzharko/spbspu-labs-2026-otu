clear; clc; close all;

base = fileparts(fileparts(fileparts(mfilename('fullpath'))));
figdir = fullfile(base, 'materials', 'figs');
resdir = fullfile(base, 'materials', 'results');
if ~exist(figdir, 'dir'), mkdir(figdir); end
if ~exist(resdir, 'dir'), mkdir(resdir); end

s = sym('s');
W1 = 1/(s + 1);
W2 = sym(2);
W3 = 5/(s^2 + 3*s + 2);
W4 = sym(4);

W = simplify((W1*W2 + W3) / (1 + W2 + W1*W2 + W3));
[N, D] = numden(W);
N = expand(N);
D = expand(D);

fid = fopen(fullfile(resdir, 'out1.txt'), 'w');
fprintf(fid, 'W|%s\n', char(W));
fprintf(fid, 'num|%s\n', strjoin(string(sym2poly(N)), ' '));
fprintf(fid, 'den|%s\n', strjoin(string(sym2poly(D)), ' '));
fclose(fid);

s_tf = tf('s');
w1 = 1/(s_tf + 1);
w2 = tf(2);
w3 = 5/(s_tf^2 + 3*s_tf + 2);
w1.InputName = 'e1'; w1.OutputName = 'y1';
w2.InputName = 'e2'; w2.OutputName = 'y2';
w3.InputName = 'e1'; w3.OutputName = 'y3';
S1 = sumblk('e1 = u - y');
S2 = sumblk('e2 = y1 - y');
S3 = sumblk('y = y2 + y3');
Sys = connect(w1, w2, w3, S1, S2, S3, 'u', 'y');
Weq = tf([2 9], [3 11 15]);

t = linspace(0, 10, 2001).';
test(1).name = 'step';
test(2).name = 'impulse';
test(3).name = 'sine';

fid = fopen(fullfile(resdir, 'out2.csv'), 'w');
fprintf(fid, 'test,max_abs_error,rms_error\n');
for k = 1:numel(test)
    switch test(k).name
        case 'step'
            [test(k).y1, ~] = step(Sys, t);
            [test(k).y2, ~] = step(Weq, t);
        case 'impulse'
            [test(k).y1, ~] = impulse(Sys, t);
            [test(k).y2, ~] = impulse(Weq, t);
        otherwise
            u = sin(t);
            test(k).y1 = lsim(Sys, u, t);
            test(k).y2 = lsim(Weq, u, t);
    end
    test(k).err = test(k).y1 - test(k).y2;
    test(k).maxerr = max(abs(test(k).err));
    test(k).rms = sqrt(mean(test(k).err.^2));
    fprintf(fid, '%s,%.6e,%.6e\n', test(k).name, test(k).maxerr, test(k).rms);
end
fclose(fid);

f = figure('Visible', 'off', 'Color', 'white', ...
           'Units', 'centimeters', 'Position', [1 1 16 18]);
set(f, 'DefaultAxesFontName', 'Liberation Serif', ...
       'DefaultTextFontName', 'Liberation Serif', ...
       'DefaultAxesFontSize', 9);
for k = 1:numel(test)
    subplot(3, 1, k);
    plot(t, test(k).y1, 'LineWidth', 1.4, 'Color', [0.05 0.25 0.60]); hold on;
    grid on; box on;
    plot(t, test(k).y2, '--', 'LineWidth', 1.0, 'Color', [0.75 0.20 0.10]);
    ylabel('y(t)');
    if k == numel(test), xlabel('t, с'); end
    title(test(k).name);
    legend({'схема', 'W'}, 'Location', 'southeast');
    legend boxoff;
    xlim([0 10]);
end
exportgraphics(f, fullfile(figdir, 'ctrl_cases.png'), 'Resolution', 300);
close(f);

