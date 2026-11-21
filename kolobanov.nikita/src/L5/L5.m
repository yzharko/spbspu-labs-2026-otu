clear;
close all;
clc;

K = 20;
T = [0.1 0.05];
M_MAX = 1.5;
OVERSHOOT_MAX = 37;
SETTLING_BAND = 0.02;

SEARCH_RESERVE = 1e-4;
PEAK_TOLERANCE = 1e-8;
GRID_SIZE = 141;
SEARCH_TIME_POINTS = 4001;
TIME_POINTS = 20001;
OBSERVATION_FACTOR = 2;
FREQUENCY_POINTS = 1200;
DECADES = 2;

DB = 20;
HALF_TURN = 180;
RIGHT_ANGLE = HALF_TURN / 2;

s = tf('s');
W = K / (s * (T(1) * s + 1) * (T(2) * s + 1));
W_closed = feedback(W, 1);
poles_initial = pole(W_closed);
t_end = ceil(OBSERVATION_FACTOR * log(SETTLING_BAND) / max(real(poles_initial)));
t_search = linspace(0, t_end, SEARCH_TIME_POINTS).';
hurwitz = [sum(T), sum(T) - K * prod(T), K * (sum(T) - K * prod(T))];

pm_min = 2 * asind(1 / (2 * M_MAX));
plant_pm = @(w) RIGHT_ANGLE - atand(w .* T(1)) - atand(w .* T(2));
w_limit = fzero(@(w) plant_pm(w) + RIGHT_ANGLE - pm_min, ...
    [0, 1 / (min(T) * tand(pm_min / 2))]);
[w_grid, phase_part] = ndgrid(w_limit * (1:GRID_SIZE) / (GRID_SIZE + 1), ...
    (0:GRID_SIZE-1) / GRID_SIZE);
phase_min = max(0, pm_min - plant_pm(w_grid));
phase_grid = phase_min + (RIGHT_ANGLE - phase_min) .* phase_part;
ratio_grid = sqrt((1 + sind(phase_grid)) ./ (1 - sind(phase_grid)));

plant_mag = @(w) K ./ (w .* hypot(1, w .* T(1)) .* hypot(1, w .* T(2)));
corrector = @(w, r) tf([r / w, 1], [1 / (r * w), 1]) / (r * plant_mag(w));
candidate = @(w, r) feedback(W * corrector(w, r), 1);
models = arrayfun(candidate, w_grid, ratio_grid, 'UniformOutput', false);
stable = cellfun(@isstable, models);
peaks = inf(size(models));
peak_gain = @(sys) getPeakGain(sys, PEAK_TOLERANCE);
peaks(stable) = cellfun(peak_gain, models(stable));
admissible = find(peaks <= M_MAX);

measure = @(sys) stepinfo(step(sys, t_search), t_search, 1, ...
    'SettlingTimeThreshold', SETTLING_BAND - SEARCH_RESERVE);
quality = cellfun(measure, models(admissible));
settling = [quality.SettlingTime];
settling(~isfinite(settling) | [quality.Overshoot] > OVERSHOOT_MAX) = Inf;
[~, best] = min(settling);
selected = admissible(best);

w_c = w_grid(selected);
ratio = ratio_grid(selected);
phase_lead = phase_grid(selected);
w_z = w_c / ratio;
w_p = w_c * ratio;
K_target = 1 / (ratio * plant_mag(w_c));
w_corners = [1 ./ T, w_z, w_p];
w = logspace(log10(min(w_corners)) - DECADES, ...
    log10(max(w_corners)) + DECADES, FREQUENCY_POINTS).';
w = unique([w; w_corners.'; w_c]);

[magnitude, phase] = bode(W, w);
L = DB * log10(squeeze(magnitude));
phi = squeeze(phase);
L_target = L + DB * log10(K_target * hypot(1, w / w_z) ./ hypot(1, w / w_p));
phi_target = phi + atand(w / w_z) - atand(w / w_p);
L_asym = DB * (log10(K ./ w) - sum(max(log10(w .* T), 0), 2));
L_target_asym = L_asym + DB * (log10(K_target) + ...
    max(log10(w / w_z), 0) - max(log10(w / w_p), 0));
L_corrector = L_target_asym - L_asym;

K_corrector = 10 ^ (L_corrector(1) / DB);
W_corrector = K_corrector * (s / w_z + 1) / (s / w_p + 1);
W_final = W_corrector * W;
W_final_closed = feedback(W_final, 1);

poles = {poles_initial, pole(W_final_closed)};
t = linspace(0, t_end, TIME_POINTS).';
h = [step(W_closed, t), step(W_final_closed, t)];
step_quality(1) = stepinfo(h(:, 1), t, dcgain(W_closed), ...
    'SettlingTimeThreshold', SETTLING_BAND);
step_quality(2) = stepinfo(h(:, 2), t, dcgain(W_final_closed), ...
    'SettlingTimeThreshold', SETTLING_BAND);
[M(1), w_peak(1)] = peak_gain(W_closed);
[M(2), w_peak(2)] = peak_gain(W_final_closed);
[gain_margin(1), phase_margin(1), w_phase(1), w_gain(1)] = margin(W);
[gain_margin(2), phase_margin(2), w_phase(2), w_gain(2)] = margin(W_final);
L_at_crossover = DB * log10(abs(freqresp(W_final, w_c)));

disp('Определители Гурвица исходной замкнутой САУ:');
disp(hurwitz);
disp('Полюса исходной и скорректированной замкнутых САУ:');
disp(poles{1});
disp(poles{2});
fprintf('Предел частоты среза %.6f рад/с, нижняя оценка запаса %.6f град.\n', ...
    w_limit, pm_min);
fprintf('Частоты среза и сопряжения: %.6f, %.6f, %.6f рад/с.\n', w_c, w_z, w_p);
fprintf('ПКУ: K = %.8f, Tz = %.8f с, Tp = %.8f с, опережение %.6f град.\n', ...
    K_corrector, 1 / w_z, 1 / w_p, phase_lead);
fprintf('ЛАЧХ скорректированной разомкнутой САУ на частоте среза: %.3g дБ.\n', ...
    L_at_crossover);
disp(table(gain_margin.', phase_margin.', w_phase.', w_gain.', M.', w_peak.', ...
    [step_quality.SettlingTime].', [step_quality.Overshoot].', ...
    'VariableNames', { ...
        'Запас по амплитуде, раз', 'Запас по фазе, град.', 'Частота фазы −180°, рад/с', ...
        'Частота среза, рад/с', 'M', 'Частота пика, рад/с', ...
        'Время установления, с', 'Перерегулирование, %'}, ...
    'RowNames', {'Исходная', 'Скорректированная'}));

open_window = figure('Name', 'Исходные и желаемые частотные характеристики', 'Color', 'w');
subplot(2, 1, 1);
semilogx(w, [L L_target]);
grid on;
ylabel('L, дБ');
yline(0, 'k:');
xline(w_c, 'k:');
legend('Исходная', 'Желаемая', 'Location', 'southwest');

subplot(2, 1, 2);
semilogx(w, [phi phi_target]);
grid on;
ylabel('Фаза, град.');
xlabel('\omega, рад/с');
xline(w_c, 'k:');

asym_window = figure('Name', 'Геометрическое вычитание ЛАЧХ', 'Color', 'w');
semilogx(w, [L_asym L_target_asym L_corrector]);
grid on;
xlabel('\omega, рад/с');
ylabel('L, дБ');
legend('Исходная', 'Желаемая', 'ПКУ', 'Location', 'southwest');
xline(w_z, 'k:', '\omega_z', 'HandleVisibility', 'off');
xline(w_p, 'k:', '\omega_p', 'HandleVisibility', 'off');

closed_window = figure('Name', 'Частотные характеристики замкнутых САУ', 'Color', 'w');
[mag_initial, phase_initial] = bode(W_closed, w);
[mag_final, phase_final] = bode(W_final_closed, w);
subplot(2, 1, 1);
semilogx(w, DB * log10([squeeze(mag_initial) squeeze(mag_final)]));
grid on;
ylabel('L, дБ');
legend('Исходная', 'Скорректированная', 'Location', 'southwest');

subplot(2, 1, 2);
semilogx(w, [squeeze(phase_initial) squeeze(phase_final)]);
grid on;
xlabel('\omega, рад/с');
ylabel('Фаза, град.');

step_window = figure('Name', 'Переходные характеристики замкнутых САУ', 'Color', 'w');
plot(t, h);
grid on;
xlabel('Время, с');
ylabel('h(t)');
xlim([0 ceil(max([step_quality.SettlingTime]))]);
legend('Исходная', 'Скорректированная', 'Location', 'northeast');
yline(1 + SETTLING_BAND, 'k:', 'HandleVisibility', 'off');
yline(1 - SETTLING_BAND, 'k:', 'HandleVisibility', 'off');
