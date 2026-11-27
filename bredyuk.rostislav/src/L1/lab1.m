% Лабораторная работа 1

clear; close all; clc;


K     = 1.7;
T     = 0.3;
T1    = 0.4;
T2    = 0.6;
T1_ap = 2*T2;

t_end   = 10;
N       = 2000;
tv      = linspace(0, t_end, N);
eps_val = 1e-3;

syms s t

W_sym = [ ...
    K/s
    K*s
    K
    K/(T*s + 1)
    K/(T2^2*s^2 + T1_ap*s + 1)
    K/(T2^2*s^2 + T1*s + 1)
    K/(T*s^2 + 1)
    K/(s*(T*s + 1))
    K*s/(T*s + 1)
    K*(T*s + 1)
    K*(T*s + 1)/s ];

names = { ...
    '1. Интегрирующее (идеальное)'; ...
    '2. Дифференцирующее (идеальное)'; ...
    '3. Усилительное (безынерционное)'; ...
    '4. Апериодическое 1-го порядка'; ...
    '5. Апериодическое 2-го порядка'; ...
    '6. Колебательное'; ...
    '7. Консервативное'; ...
    '8. Инерционное интегрирующее'; ...
    '9. Реальное дифференцирующее'; ...
    '10. Форсирующее 1-го порядка'; ...
    '11. Изодромное' };

%% --------------- 3. Аналитический расчёт h(t) и w(t) ---------------
U_step = 1/s;      % изображение единичного скачка 1(t)
U_imp  = sym(1);   % изображение единичного импульса delta(t)

h_ex = cell(11, 1);
w_ex = cell(11, 1);
for i = 1:11
    h_ex{i} = simplify(ilaplace(W_sym(i)*U_step, s, t));
    w_ex{i} = simplify(ilaplace(W_sym(i)*U_imp,  s, t));
end

fprintf('Аналитические характеристики звеньев:\n');
for i = 1:11
    fprintf('%s\n  h(t) = %s\n  w(t) = %s\n', names{i}, char(h_ex{i}), char(w_ex{i}));
end

W_tf = { ...
    tf(K, [1 0])            % 1:  K/s
    tf(K*[1 0], [eps_val 1]) % 2:  K*s/(eps_val*s+1)
    tf(K, 1)                % 3:  K
    tf(K, [T 1])            % 4:  K/(T*s+1)
    tf(K, [T2^2 T1_ap 1])   % 5:  K/(T2^2*s^2+T1_ap*s+1)
    tf(K, [T2^2 T1 1])      % 6:  K/(T2^2*s^2+T1*s+1)
    tf(K, [T 0 1])          % 7:  K/(T*s^2+1)
    tf(K, [T 1 0])          % 8:  K/(T*s^2+s)
    tf(K*[1 0], [T 1])      % 9:  K*s/(T*s+1)
    tf(K*[T 1], [eps_val 1]) % 10: K*(T*s+1)/(eps_val*s+1)
    tf(K*[T 1], [1 0]) };   % 11: (K*T*s+K)/s


fprintf('\nПроверка аналитики численным моделированием:\n');
idx_h = [1 3 4 5 6 7 8 9 11];   % звенья без дельта-составляющей в h(t)
for i = idx_h
    y_num = step(W_tf{i}, tv);
    y_an  = eval_reg(h_ex{i}, tv);
    dh = max(abs(y_num(:) - y_an(:)));
    fprintf('звено %2d:  max|h_num - h_an| = %.2e\n', i, dh);
end
idx_w = [1 4 5 6 7 8];           % звенья без дельта-составляющей в w(t)
for i = idx_w
    y_num = impulse(W_tf{i}, tv);
    y_an  = eval_reg(w_ex{i}, tv);
    dw = max(abs(y_num(:) - y_an(:)));
    fprintf('звено %2d:  max|w_num - w_an| = %.2e\n', i, dw);
end

plot_grid(h_ex, 'h(t)', 'Переходные характеристики h(t)', tv, names);
plot_grid(w_ex, 'w(t)', 'Импульсные характеристики w(t) (регулярные составляющие)', tv, names);

mult = [1/5, 1/2, 1, 2, 5];
lbl4 = {'T/5 = 0,06 с', 'T/2 = 0,15 с', 'T = 0,3 с', '2T = 0,6 с', '5T = 1,5 с'};
sty  = {'-o', '-s', '-d', '-^', '-v'};

figure('Color', 'w'); hold on; grid on;
for k = 1:numel(mult)
    Tk = T*mult(k);
    hk = K*(1 - exp(-tv/Tk));
    plot(tv, hk, sty{k}, 'LineWidth', 1.2, 'MarkerIndices', 1:100:N, ...
         'DisplayName', lbl4{k});
end
xlabel('t, с'); ylabel('h(t)');
title('Звено №4: влияние постоянной времени T');
legend('show', 'Location', 'southeast');

T1_set = [0.4, 0.2, 0.1, 0.05];
figure('Color', 'w'); hold on; grid on;
for k = 1:numel(T1_set)
    sigma = T1_set(k)/(2*T2^2);          % коэффициент затухания
    wd    = sqrt(1/T2^2 - sigma^2);      % частота затухающих колебаний
    hk    = K*(1 - exp(-sigma*tv).*(cos(wd*tv) + sigma/wd*sin(wd*tv)));
    xi    = T1_set(k)/(2*T2);            % степень затухания (демпфирования)
    overs = exp(-pi*xi/sqrt(1 - xi^2));  % перерегулирование (относительное)
    plot(tv, hk, 'LineWidth', 1.2, ...
         'DisplayName', sprintf('T_1 = %.2f с (\\xi = %.3f, \\sigma = %.0f%%)', ...
                                T1_set(k), xi, 100*overs));
end
xlabel('t, с'); ylabel('h(t)');
title(sprintf('Звено №6: влияние T_1 при фиксированном T_2 = %.1f с', T2));
legend('show');

function y = eval_reg(expr, tv)
    syms t
    expr = subs(expr, dirac(1, t), sym(0));   % производная дельта-функции
    expr = subs(expr, dirac(t),    sym(0));   % дельта-функция
    f    = matlabFunction(expr, 'vars', t);
    y    = f(tv);
    y    = reshape(y, size(tv));
end

function plot_grid(exprs, ylab, ttl, tv, names)
    figure('Color', 'w');
    tiledlayout(4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    for i = 1:11
        nexttile; hold on; box on;
        plot(tv, eval_reg(exprs{i}, tv), 'LineWidth', 1.2, 'Color', [0 0.45 0.74]);
        grid on;
        xlabel('t, с'); ylabel(ylab);
        title(names{i}, 'Interpreter', 'none', 'FontSize', 9);
    end
    sgtitle(ttl);
end
