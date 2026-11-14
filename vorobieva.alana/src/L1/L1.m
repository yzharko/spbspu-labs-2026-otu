clear; clc; close all;

set(0, 'DefaultFigurePosition', [50 50 1400 900]);

K = 3.8;
T = 1.1;
T_1 = 1.5;
T_2 = 1.6;
T_1a = 2 * T_2;

syms s t k

names = { ...
  'Интегрирующее', 'Дифференцирующее', 'Усилительное', ...
  'Апериодическое 1-го порядка', 'Апериодическое 2-го порядка', ...
  'Колебательное', 'Консервативное', 'Инерционное интегрирующее', ...
  'Реальное дифференцирующее', 'Форсирующее', 'Изодромное' };

n = numel(names);

W = { ...
  k / s, ...
  k * s, ...
  k, ...
  k / (T * s + 1), ...
  k / (T_2 ^ 2 * s ^ 2 + T_1a * s + 1), ...
  k / (T_2 ^ 2 * s ^ 2 + T_1 * s + 1), ...
  k / (T * s ^ 2 + 1), ...
  k / (s * (T * s + 1)), ...
  k * s / (T * s + 1), ...
  k * (T * s + 1), ...
  k * (T * s + 1) / s };

U_step = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);

h = cell(1, n);
w = cell(1, n);

for i = 1:n
  Y_step = W{ i } * U_step;
  Y_impulse = W{ i } * U_impulse;

  h{ i } = simplify(ilaplace(Y_step, s, t));
  w{ i } = simplify(ilaplace(Y_impulse, s, t));

  fprintf('%d. %s\n   W(s) = %s\n   h(t) = %s\n   w(t) = %s\n', i, names{ i }, char(W{ i }), char(h{ i }), char(w{ i }));
end

s_tf = tf('s');

W_tf = { K / s_tf, K * s_tf, tf(K), K / (T * s_tf + 1), ...
  K / (T_2 ^ 2 * s_tf ^ 2 + T_1a * s_tf + 1), ...
  K / (T_2 ^ 2 * s_tf ^ 2 + T_1 * s_tf + 1), K / (T * s_tf ^ 2 + 1), ...
  K / (s_tf * (T * s_tf + 1)), K * s_tf / (T * s_tf + 1), ...
  K * (T * s_tf + 1), K * (T * s_tf + 1) / s_tf };

time = 0:0.01:15;
f = { h, w };
y_names = { 'h(t)', 'w(t)' };
fig_names = { 'Переходные характеристики', 'Импульсные характеристики' };

for m = 1:2
  figure('Name', fig_names{ m });
  
  for i = 1:n
    subplot(3, 4, i); 
    hold on; 
    grid on;
    y = subs(f{ m }{ i }, [dirac(t) dirac(1, t)], [0 0]);
    fplot(subs(y, k, K), [0 15], 'LineWidth', 2);

    if isproper(W_tf{ i }) && order(W_tf{ i }) > 0
      if m == 1
        plot(time, step(W_tf{ i }, time), '--');
      else
        plot(time, impulse(W_tf{ i }, time), '--');
      end
    end

    if contains(char(f{ m }{ i }), 'dirac')
      delta_part = simplify(subs(f{ m }{ i } - y, k, K));
      text(0.05, 0.9, ['+ ' char(vpa(delta_part, 3))], ...
                'Units', 'normalized', 'Color', 'r');
    end
    title(sprintf('%d. %s', i, names{ i }), 'FontSize', 9);
    xlabel('Время t, с');
    ylabel(y_names{ m });
  end
  sgtitle([fig_names{ m } ' типовых звеньев']);
end

t_T  = 0:0.01:30;
t_T1 = 0:0.01:80;
K_values = K * [0.5 1 2];
sys_K = { @(Ki) Ki / (T * s_tf + 1), ...
          @(Ki) Ki / (T_2 ^ 2 * s_tf ^ 2 + T_1 * s_tf + 1) };
sys_numbers = [4 6];

figure('Name', 'Влияние K');
for j = 1:2
  subplot(1, 2, j); 
  hold on; 
  grid on;
  
  for Ki = K_values
    plot(time, step(sys_K{ j }(Ki), time), 'DisplayName', sprintf('K = %g', Ki));
  end
  
  title(sprintf('%d. %s', sys_numbers(j), names{ sys_numbers(j) }));
  xlabel('Время t, с');
  ylabel('h(t)');
  legend('show', 'Location', 'southeast');
end

T_values = T * [0.2 0.5 1 2 5];

figure('Name', 'Влияние T, апериодическое звено');
hold on; 
grid on;
fprintf('\nАпериодическое звено (K = %g)\n', K);
fprintf('   T      t95%%, с   (время выхода на 95%% от K)\n');

for Ti = T_values
  y = step(K / (Ti * s_tf + 1), t_T);
  plot(t_T, y, 'DisplayName', sprintf('T = %.2g с', Ti));
  idx = find(y >= 0.95 * K, 1);
  fprintf('%6.2f  %7.2f\n', Ti, t_T(idx));
end

title('Апериодическое звено при разных T');
xlabel('Время t, с');
ylabel('h(t)');
xlim([0 30]);
legend('show', 'Location', 'southeast');
hold off;

T1_values = T_1 * [1 0.6 0.4 0.2 0];

figure('Name', 'Влияние T1, колебательное звено');
hold on; 
grid on;
fprintf('\nКолебательное звено (K = %g, T_2 = %g)\n', K, T_2);
fprintf('   T1     zeta   sigma_теор, %%   sigma_график, %%   t5%%, с\n');

for T1i = T1_values
  y = step(K / (T_2 ^ 2 * s_tf ^ 2 + T1i * s_tf + 1), t_T1);
  plot(t_T1, y, 'DisplayName', sprintf('T1 = %.2g с', T1i));
  overshoot = (max(y) - K) / K * 100;
  idx = find(abs(y - K) > 0.05 * K, 1, 'last');
  if idx == numel(t_T1)
    t_settle = NaN;
  else
    t_settle = t_T1(idx + 1);
  end
  zeta = T1i / (2 * T_2); 
  sigma_theory = 100 * exp(-pi * zeta / sqrt(1 - zeta ^ 2));
  fprintf('%6.2f  %6.3f  %10.2f  %13.2f  %10.2f\n', T1i, zeta, sigma_theory, overshoot, t_settle);
end

title('Колебательное звено при уменьшении T1');
xlabel('Время t, с');
ylabel('h(t)');
xlim([0 80]);
legend('show', 'Location', 'eastoutside');
hold off;

figs = findobj('Type', 'figure');

for j = 1:numel(figs)
  exportgraphics(figs(j), sprintf('figure_%d.png', figs(j).Number), 'Resolution', 200);
end
