clear;
clc;
close all;

K0 = 3.2;
T0 = 0.9;
T1 = 1.0;
T2 = 1.4;
T1a = 3.0;

syms s t
syms K T positive

names = {
  'Интегрирующее'
  'Идеальное дифференцирующее'
  'Усилительное'
  'Апериодическое 1-го порядка'
  'Апериодическое 2-го порядка'
  'Колебательное'
  'Консервативное'
  'Инерционное интегрирующее'
  'Реальное дифференцирующее'
  'Форсирующее 1-го порядка'
  'Изодромное'
};

W = [
  K / s;
  K * s;
  K;
  K / (T * s + 1);
  K / (sym(T2)^2 * s^2 + T1a * s + 1);
  K / (sym(T2)^2 * s^2 + T1 * s + 1);
  K / (T * s^2 + 1);
  K / (s * (T * s + 1));
  K * s / (T * s + 1);
  K * (T * s + 1);
  K * (T * s + 1) / s;
];

N = numel(W);

U_step = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);

h = sym(zeros(N, 1));
w = sym(zeros(N, 1));

for i = 1:N
  Y_step = simplify(W(i) * U_step);
  Y_impulse = simplify(W(i) * U_impulse);

  h(i) = simplify(ilaplace(Y_step, s, t));
  w(i) = simplify(ilaplace(Y_impulse, s, t));

  fprintf('%d. %s\n  Y_step = %s\n  Y_impulse = %s\n  h(t) = %s\n  w(t) = %s\n', ...
    i, names{i}, char(Y_step), char(Y_impulse), char(h(i)), char(w(i)));
end

responses = subs([h w], [K T], [K0 T0]);
regular = subs(subs(responses, dirac(1, t), 0), dirac(t), 0);

delta_coef = double(subs(subs(responses, dirac(1, t), 0), dirac(t), 1) - regular);
delta_prime_coef = double(subs(subs(responses, dirac(t), 0), dirac(1, t), 1) - regular);

regular = subs(regular, heaviside(t), 1);

W_numeric = subs(W, [K T], [K0 T0]);
time = 0:0.01:20;

figures = gobjects(5, 1);
labels = {'h(t)', 'w(t)'};
titles = {'Переходные характеристики h(t)', 'Импульсные характеристики w(t)'};
deviation = zeros(N, 2);

for k = 1:2
  figures(k) = figure('Name', titles{k}, 'Position', [50 50 1100 1080]);
  tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

  for i = 1:N
    nexttile;
    hold on;
    grid on;

    y = double(subs(regular(i, k), t, time)) + zeros(size(time));

    [num, den] = numden(W_numeric(i));
    num = sym2poly(num);
    den = sym2poly(den);

    if k == 1
      den = conv(den, [1 0]);
    end

    if numel(num) >= numel(den)
      [q, r] = deconv(num, den);
    else
      q = 0;
      r = num;
    end

    q = [zeros(1, 2 - numel(q)), q];

    assert(max(abs(q - [delta_prime_coef(i, k), delta_coef(i, k)])) < 1e-9, ...
      'Не совпали коэффициенты импульсов у звена %d.', i);

    if any(abs(r) > 1e-12)
      y_sim = reshape(impulse(tf(r, den), time), 1, []);
    else
      y_sim = zeros(size(time));
    end

    deviation(i, k) = max(abs(y - y_sim));

    assert(deviation(i, k) < 1e-7 * max(1, max(abs(y))), ...
      'Не совпали обычные части отклика у звена %d.', i);

    plot(time, y, 'LineWidth', 2);
    plot(time, y_sim, '--', 'LineWidth', 1);

    note = '';

    if delta_coef(i, k) ~= 0
      note = sprintf('%+.4g delta(t)', delta_coef(i, k));
    end

    if delta_prime_coef(i, k) ~= 0
      note = [note, sprintf(' %+.4g delta''(t)', delta_prime_coef(i, k))];
    end

    text(0.02, 0.95, note, 'Units', 'normalized', ...
      'VerticalAlignment', 'top', 'Color', 'r', 'Interpreter', 'none');

    title(sprintf('%d. %s', i, names{i}));
    xlabel('Время t, с');
    ylabel(labels{k});
    xlim([0 time(end)]);

    if i == 1
      legend('аналитический расчёт', 'моделирование', 'Location', 'best');
    end
  end

  sgtitle(sprintf('%s при K=%g, T=%g, T2=%g; T1(5)=%g, T1(6)=%g', ...
    titles{k}, K0, T0, T2, T1a, T1));
end

disp(table((1:N)', deviation(:, 1), deviation(:, 2), ...
  'VariableNames', {'Link', 'Error_h', 'Error_w'}));

s_tf = tf('s');

figures(3) = figure('Name', 'Влияние K', 'Position', [50 50 1100 450]);
tiledlayout(1, 2);

ax4 = nexttile;
hold on;
grid on;

ax6 = nexttile;
hold on;
grid on;

for Ki = K0 * [0.5 1 2]
  plot(ax4, time, step(Ki / (T0 * s_tf + 1), time), ...
    'DisplayName', sprintf('K = %g', Ki));

  plot(ax6, time, step(Ki / (T2^2 * s_tf^2 + T1 * s_tf + 1), time), ...
    'DisplayName', sprintf('K = %g', Ki));
end

title(ax4, '4. Апериодическое 1-го порядка');
title(ax6, '6. Колебательное');

for ax = [ax4 ax6]
  xlabel(ax, 'Время t, с');
  ylabel(ax, 'h(t)');
  legend(ax, 'Location', 'best');
end

time_aperiodic = 0:0.01:40;

figures(4) = figure('Name', 'Влияние T', 'Position', [50 50 1100 450]);
tiledlayout(1, 2);

ax_h = nexttile;
hold on;
grid on;

ax_w = nexttile;
hold on;
grid on;

fprintf('\nЗвено 4: T, время установления (5%%), w(0+)\n');

for Ti = T0 * [0.2 0.5 1 2 5]
  model = K0 / (Ti * s_tf + 1);
  y_h = step(model, time_aperiodic);
  y_w = impulse(model, time_aperiodic);

  plot(ax_h, time_aperiodic, y_h, 'DisplayName', sprintf('T = %g', Ti));
  plot(ax_w, time_aperiodic, y_w, 'DisplayName', sprintf('T = %g', Ti));

  fprintf('%6g %10.3f %10.3f\n', Ti, -Ti * log(0.05), y_w(1));
end

sgtitle(sprintf('Звено 4 при K=%g и различных T', K0));

for ax = [ax_h ax_w]
  xlabel(ax, 'Время t, с');
  legend(ax, 'Location', 'best');
end

ylabel(ax_h, 'h(t)');
ylabel(ax_w, 'w(t)');

time_oscillatory = 0:0.01:80;

figures(5) = figure('Name', 'Влияние T1', 'Position', [50 50 1100 700]);
tiledlayout(2, 1);

ax_h = nexttile;
hold on;
grid on;

ax_w = nexttile;
hold on;
grid on;

fprintf('\nЗвено 6: T1, перерегулирование (%%), время установления (5%%)\n');

for T1i = T1 * [1 0.75 0.5 0.25 0]
  model = K0 / (T2^2 * s_tf^2 + T1i * s_tf + 1);
  y_h = step(model, time_oscillatory);
  y_w = impulse(model, time_oscillatory);

  plot(ax_h, time_oscillatory, y_h, 'DisplayName', sprintf('T1 = %g', T1i));
  plot(ax_w, time_oscillatory, y_w, 'DisplayName', sprintf('T1 = %g', T1i));

  if T1i == 0
    overshoot = NaN;
    settling_time = Inf;
  else
    zeta = T1i / (2 * T2);
    overshoot = 100 * exp(-pi * zeta / sqrt(1 - zeta^2));
    last = find(abs(y_h - K0) > 0.05 * K0, 1, 'last');
    settling_time = time_oscillatory(last + 1);
  end

  fprintf('%6g %10.2f %10.2f\n', T1i, overshoot, settling_time);
end

sgtitle(sprintf('Звено 6 при K=%g, T2=%g и уменьшении T1', K0, T2));

for ax = [ax_h ax_w]
  xlabel(ax, 'Время t, с');
  xlim(ax, [0 80]);
  legend(ax, 'Location', 'eastoutside');
end

ylabel(ax_h, 'h(t)');
ylabel(ax_w, 'w(t)');

folder = fullfile(fileparts(mfilename('fullpath')), 'results_L1_clean');

if ~exist(folder, 'dir')
  mkdir(folder);
end

files = {
  'step_overview'
  'impulse_overview'
  'K_variation'
  'T_variation'
  'T1_variation'
};

for i = 1:numel(figures)
  exportgraphics(figures(i), fullfile(folder, [files{i} '.png']), 'Resolution', 150);
end
