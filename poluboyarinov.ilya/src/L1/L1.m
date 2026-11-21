clear;
clc;
close all;

K = 1.8;
T = 0.6;
T1 = 0.9;
T2 = 1.2;
T11 = 3.0; %чтобы w5 и w6 не были одинаковыми

s = tf('s');


W1 = K / s;

W2 = K * s;

W3 = K;

W4 = K / (T*s + 1);

W5 = K / (T2^2*s^2 + T11*s + 1);

W6 = K / (T2^2*s^2 + T1*s + 1);

W7 = K / (T*s^2 + 1);

W8 = K / (s*(T*s + 1));

W9 = K*s / (T*s + 1);

W10 = K * (T*s + 1);

W11 = K * (T*s + 1) / s;


disp('Исходные параметры:');
fprintf('K = %.2f, T = %.2f, T1 = %.2f, T2 = %.2f\n', K, T, T1, T2);

disp('Передаточные функции:');

for i = 1:11
    fprintf('\nЗвено №%d:\n', i);
    eval(['W = W' num2str(i) ';']);
    disp(W);
end

if T11 >= 2*T2
    disp('Звено №5: апериодическое второго порядка');
else
    disp('Звено №5: колебательное');
end

if T1 >= 2*T2
    disp('Звено №6: апериодическое второго порядка');
else
    disp('Звено №6: колебательное');
end


syms p tau k real

Ws = {
    k/p, ...
    k*p, ...
    k, ...
    k/(T*p + 1), ...
    k/(T2^2*p^2 + T11*p + 1), ...
    k/(T2^2*p^2 + T1*p + 1), ...
    k/(T*p^2 + 1), ...
    k/(p*(T*p + 1)), ...
    k*p/(T*p + 1), ...
    k*(T*p + 1), ...
    k*(T*p + 1)/p
};

U_step = 1/p;
U_impulse = 1;

Y_step = cell(11,1);
Y_impulse = cell(11,1);
h = cell(11,1);
w = cell(11,1);

for i = 1:11

    Y_step{i} = simplify(Ws{i} * U_step);

    Y_impulse{i} = simplify(Ws{i} * U_impulse);

    h{i} = simplify(ilaplace(Y_step{i}, p, tau));
    w{i} = simplify(ilaplace(Y_impulse{i}, p, tau));

    fprintf('\n====================================\n');
    fprintf('ЗВЕНО №%d\n', i);
    fprintf('====================================\n');

    fprintf('W(s) = \n');
    disp(Ws{i});

    fprintf('Y_step(s) = \n');
    disp(Y_step{i});

    fprintf('h(t) = \n');
    disp(h{i});

    fprintf('Y_impulse(s) = \n');
    disp(Y_impulse{i});

    fprintf('w(t) = \n');
    disp(w{i});

end


K_values = [0.9 1.8 3.6];

time = 0:0.01:10;

names = {'Интегрирующее', ...
    'Идеальное дифференцирующее', ...
    'Усилительное', ...
    'Апериодическое 1-го порядка', ...
    'Апериодическое 2-го порядка', ...
    'Колебательное', ...
    'Консервативное', ...
    'Инерционное интегрирующее', ...
    'Реальное дифференцирующее', ...
    'Форсирующее', ...
    'Изодромное'};

notes_h = cell(11,1);
notes_w = cell(11,1);

notes_h{2} = 'Содержит K*delta(t)';
notes_h{10} = 'Содержит K*T*delta(t)';

notes_w{2} = 'Содержит K*delta''(t)';
notes_w{3} = 'Содержит K*delta(t)';
notes_w{9} = 'Содержит импульсную составляющую';
notes_w{10} = 'Содержит delta(t) и delta''(t)';
notes_w{11} = 'Содержит K*T*delta(t)';

figure('Name','Переходные характеристики', ...
    'Color','w','Position',[50 50 1450 900]);

for i = 1:11

    subplot(4,3,i);
    hold on;
    grid on;

    for K_plot = K_values

        expr = subs(h{i},k,K_plot);
        expr = subs(expr,dirac(1,tau),0);
        expr = subs(expr,dirac(tau),0);
        expr = subs(expr,heaviside(tau),1);

        f = matlabFunction(expr,'Vars',tau);
        y = f(time);

        if isscalar(y)
            y = repmat(y,size(time));
        end

        plot(time,y,'LineWidth',1.3, ...
            'DisplayName',sprintf('K = %.1f',K_plot));
    end

    title(sprintf('%d. %s',i,names{i}));
    xlabel('Время, с');
    ylabel('h(t)');

    if ~isempty(notes_h{i})
        text(0.05,0.90,notes_h{i}, ...
            'Units','normalized','FontSize',8);
    end

    if i == 1
        legend('show','Location','best');
    end
end

sgtitle('Переходные характеристики динамических звеньев');

figure('Name','Импульсные характеристики', ...
    'Color','w','Position',[50 50 1450 900]);

for i = 1:11
    subplot(4, 3, i);
    hold on;
    grid on;

    for K_plot = K_values
        expr = subs(w{i}, k, K_plot);

        expr = subs(expr, dirac(1, tau), 0);
        expr = subs(expr, dirac(tau), 0);

        expr = subs(expr, heaviside(tau), 1);

        f = matlabFunction(expr, 'Vars', tau);

        y = f(time);

        if isscalar(y)
            y = repmat(y, size(time));
        end

        plot(time, y, 'LineWidth', 1.3, ...
            'DisplayName', sprintf('K = %.1f', K_plot));
    end

    title(sprintf('Звено №%d', i));
    xlabel('Время t');
    ylabel('w(t)');

    if i == 1
        legend('show', 'Location', 'best');
    end

    hold off;
end

sgtitle('Импульсные характеристики динамических звеньев');

Wtf = {W1, W2, W3, W4, W5, W6, W7, W8, W9, W10, W11};

time_check = linspace(0, 10, 1000);
step_check = [1 3 4 5 6 7 8 9 11];

fprintf('\nПРОВЕРКА ПЕРЕХОДНЫХ ХАРАКТЕРИСТИК\n');

for i = step_check

    if i == 3
        y_analytic = K * ones(size(time_check));
        y_matlab = K * ones(size(time_check));

        error_value = max(abs(y_matlab - y_analytic));

        fprintf('Звено №3: ошибка = %.10e\n', error_value);
        continue;
    end

    expr = subs(h{i}, k, K);
    expr = subs(expr, dirac(1, tau), 0);
    expr = subs(expr, dirac(tau), 0);
    expr = subs(expr, heaviside(tau), 1);

    f = matlabFunction(expr, 'Vars', tau);
    y_analytic = f(time_check);
    y_analytic = y_analytic(:);

    y_matlab = step(Wtf{i}, time_check);
    y_matlab = squeeze(y_matlab);
    y_matlab = y_matlab(:);

    if numel(y_matlab) ~= numel(y_analytic)
        fprintf('Звено №%d: несовпадение размеров массивов\n', i);
        continue;
    end

    error_value = max(abs(y_matlab - y_analytic));

    fprintf('Звено №%d: ошибка = %.10e\n', i, error_value);

end

impulse_check = [1 4 5 6 7 8];

fprintf('\nПРОВЕРКА ИМПУЛЬСНЫХ ХАРАКТЕРИСТИК\n');

for i = impulse_check

    expr = subs(w{i}, k, K);
    expr = subs(expr, dirac(1, tau), 0);
    expr = subs(expr, dirac(tau), 0);
    expr = subs(expr, heaviside(tau), 1);

    f = matlabFunction(expr, 'Vars', tau);
    y_analytic = f(time_check);

    if isscalar(y_analytic)
        y_analytic = repmat(y_analytic, size(time_check));
    end

    y_matlab = impulse(Wtf{i}, time_check);

    error_value = max(abs(y_matlab(:) - y_analytic(:)));

    fprintf('Звено №%d: ошибка = %.10e\n', i, error_value);

end

K_fixed = K;
T_values = [0.12 0.3 0.6 1.2 3.0];

time_T = linspace(0, 15, 1500);

figure('Name', 'Влияние T на переходную характеристику');

hold on;
grid on;

for T_plot = T_values

    h_T = K_fixed * (1 - exp(-time_T / T_plot));

    plot(time_T, h_T, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('T = %.2f', T_plot));

end

title('Влияние постоянной времени T на h(t)');
xlabel('Время t');
ylabel('h(t)');
legend('show', 'Location', 'best');

hold off;


figure('Name', 'Влияние T на импульсную характеристику');

hold on;
grid on;

for T_plot = T_values

    w_T = (K_fixed / T_plot) * exp(-time_T / T_plot);

    plot(time_T, w_T, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('T = %.2f', T_plot));

end

title('Влияние постоянной времени T на w(t)');
xlabel('Время t');
ylabel('w(t)');
legend('show', 'Location', 'best');

hold off;