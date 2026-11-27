clear;
clc;
close all;

SLOW_POLE_RATE = 0.2;
FAST_POLE_RATE = 4;
MAX_OVERSHOOT = 8;
MAX_SETTLING_TIME = 5;
SETTLING_BAND = 0.02;
ZN_FACTORS = [0.6 0.5 0.125];
PLOT_SIZE = [100 100 900 440];
FONT_SIZE = 14;
LINE_WIDTH = 1.6;

s = tf('s');
Go = 1 / ((s + SLOW_POLE_RATE)^2 * (s + FAST_POLE_RATE));
[~, den] = tfdata(Go, 'v');
G0 = feedback(Go, 1);
initial_poles = pole(G0);
[Kcr, phase_margin_initial, wcr] = margin(Go);
Tcr = 2 * pi / wcr;

zn_settings = ZN_FACTORS .* [Kcr Tcr Tcr];
Czn = pid(zn_settings(1), zn_settings(1) / zn_settings(2), zn_settings(1) * zn_settings(3));
Cauto = pidtune(Go, 'PID');
Lzn = Czn * Go;
Lauto = Cauto * Go;
Gzn = feedback(Lzn, 1);
Gauto = feedback(Lauto, 1);

models = [G0; Gzn; Gauto];
[response, time] = step(models);
response = real(response);
steady_output = dcgain(models);
quality = stepinfo(response, time, steady_output, 'SettlingTimeThreshold', SETTLING_BAND);
settling_time = [quality.SettlingTime].';
overshoot = [quality.Overshoot].';
steady_error = 1 - steady_output;
stable = [isstable(G0); isstable(Gzn); isstable(Gauto)];
[gain_margins, phase_margins] = margin(stack(3, Lzn, Lauto));
gain_margin_db = 20 * log10([Kcr; gain_margins(:)]);
phase_margin = [phase_margin_initial; phase_margins(:)];
meets_time_limit = settling_time <= MAX_SETTLING_TIME;
meets_overshoot_limit = overshoot <= MAX_OVERSHOOT;

disp(initial_poles);
disp(table(Kcr, wcr, Tcr, zn_settings(2), zn_settings(3), ...
    'VariableNames', {'Kcr', 'wcr', 'Tcr', 'Ti', 'Td'}));
disp(table([Czn.Kp; Cauto.Kp], [Czn.Ki; Cauto.Ki], [Czn.Kd; Cauto.Kd], ...
    'VariableNames', {'Kp', 'Ki', 'Kd'}, 'RowNames', {'ZieglerNichols', 'pidtune'}));
disp(table(stable, settling_time, overshoot, steady_error, gain_margin_db, ...
    phase_margin, meets_time_limit, meets_overshoot_limit, ...
    'RowNames', {'Initial', 'ZieglerNichols', 'pidtune'}));

comparison_figure = figure('Name', 'Переходные характеристики', 'Color', 'w', 'Position', PLOT_SIZE, ...
    'DefaultAxesFontSize', FONT_SIZE, 'DefaultLineLineWidth', LINE_WIDTH);
plot(time, response(:, 2), '-', time, response(:, 3), '--', time, response(:, 1), 'k:');
grid on;
xlabel('Время t, с');
ylabel('Выход y(t), отн. ед.');
legend('Зиглер-Николс', 'Автонастройка pidtune', 'Без регулятора', ...
    'Location', 'southoutside', 'Orientation', 'horizontal');
