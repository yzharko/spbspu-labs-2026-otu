% Аналитическая и численная проверка варианта 8.
% Модель Simulink с тремя источниками: scheme_L2.slx.

clear; clc; close all;

here = fileparts(mfilename('fullpath'));
base = fileparts(fileparts(here));
resdir = fullfile(base, 'materials', 'results');
if ~exist(resdir, 'dir'), mkdir(resdir); end

lab = struct();
lab.variant = 8;
lab.W1 = 1 / (sym('s') + 1);
lab.W2 = sym(2);
lab.W3 = 5 / (sym('s')^2 + 3*sym('s') + 2);
lab.W4 = sym(4);

s = sym('s');
W = simplify(lab.W3 * (1 + lab.W1*lab.W2) / ...
             (1 + lab.W3 * (2 + lab.W2 + lab.W1*lab.W2)));
[N, D] = numden(W);
N = expand(N);
D = expand(D);
fprintf('Вариант %d\n', lab.variant);
fprintf('W(s) = %s\n', char(W));
fprintf('num = %s\n', char(N));
fprintf('den = %s\n', char(D));

fid = fopen(fullfile(resdir, 'analytic.txt'), 'w');
fprintf(fid, 'W|%s\n', char(W));
fprintf(fid, 'num|%s\n', strjoin(string(sym2poly(N)), ' '));
fprintf(fid, 'den|%s\n', strjoin(string(sym2poly(D)), ' '));
fclose(fid);

s_tf = tf('s');
w1 = 1 / (s_tf + 1);
w2 = tf(2);
w3 = 5 / (s_tf^2 + 3*s_tf + 2);
w1.InputName = 'e1'; w1.OutputName = 'y1';
w2.InputName = 'e2'; w2.OutputName = 'y2';
w3.InputName = 'e3'; w3.OutputName = 'y';
S1 = sumblk('e1 = u - y');
S2 = sumblk('e2 = y1 - y');
S3 = sumblk('e3 = y2 + e1 - y');
Sys = connect(w1, w2, w3, S1, S2, S3, 'u', 'y');
Weq = tf([5 15], [1 4 25 32]);

fprintf('Степень исходной модели: %d\n', numel(pole(minreal(Sys))));
fprintf('Степень эквивалентной модели: %d\n', numel(pole(Weq)));

t = linspace(0, 10, 2001).';
cases = {'step', 'impulse', 'sine'};
fid = fopen(fullfile(resdir, 'metrics_ctrl.csv'), 'w');
fprintf(fid, 'test,max_abs_error,rms_error\n');
for k = 1:numel(cases)
    switch cases{k}
        case 'step'
            [y1, ~] = step(Sys, t);
            [y2, ~] = step(Weq, t);
        case 'impulse'
            [y1, ~] = impulse(Sys, t);
            [y2, ~] = impulse(Weq, t);
        otherwise
            u = sin(t);
            y1 = lsim(Sys, u, t);
            y2 = lsim(Weq, u, t);
    end
    e = y1 - y2;
    fprintf(fid, '%s,%.6e,%.6e\n', cases{k}, max(abs(e)), sqrt(mean(e.^2)));
end
fclose(fid);

fprintf('\nГотово: %s\n', resdir);
