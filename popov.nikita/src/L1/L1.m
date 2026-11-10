clear; clc; close all;

K = 2.8; T = 0.6; T1 = 0.6; T2 = 0.5;
K_values = [K/2, K, 2*K];
syms s t
U_step = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);
fprintf('U_step(s) = %s\nU_impulse(s) = %s\n', char(U_step), char(U_impulse));

names = {
    'Идеальное интегрирующее';
    'Идеальное дифференцирующее';
    'Усилительное';
    'Апериодическое 1-го порядка';
    'Апериодическое 2-го порядка';
    'Колебательное';
    'Консервативное';
    'Инерционное интегрирующее';
    'Реальное дифференцирующее';
    'Форсирующее 1-го порядка';
    'Изодромное'
};

W = {
    K/s;
    K*s;
    sym(K);
    K/(T*s+1);
    K/(T2^2*s^2+T1*s+1);
    K/(T2^2*s^2+T1*s+1);
    K/(T*s^2+1);
    K/(s*(T*s+1));
    K*s/(T*s+1);
    K*(T*s+1);
    K*(T*s+1)/s
};

if T1 < 2*T2
    warning('Звено №5: T1 < 2*T2; при исходных данных звенья №5 и №6 совпадают.');
end

Y_step = cell(1,11); Y_impulse = cell(1,11);
h = cell(1,11); w = cell(1,11);
for i = 1:11
    Y_step{i} = simplify(W{i}*U_step);
    Y_impulse{i} = simplify(W{i}*U_impulse);
    h{i} = simplify(ilaplace(Y_step{i}, s, t));
    w{i} = simplify(ilaplace(Y_impulse{i}, s, t));

    fprintf('\n========== Звено №%d: %s ==========\n', i, names{i});
    fprintf('W(s)         = %s\n', char(W{i}));
    fprintf('Y_step(s)    = %s\n', char(Y_step{i}));
    fprintf('h(t)         = %s\n', char(h{i}));
    fprintf('Y_impulse(s) = %s\n', char(Y_impulse{i}));
    fprintf('w(t)         = %s\n', char(w{i}));
end

time = 0:0.01:15;
H = cell(1,11); I = cell(1,11);
for i = 1:11
    H{i} = regular_values(h{i}, t, time);
    I{i} = regular_values(w{i}, t, time);
end

noteH = cell(1,11); noteI = cell(1,11);
K_symbolic = sym('K');
for i = 1:11
    h_regular = subs(subs(h{i}, dirac(1,t), sym(0)), dirac(t), sym(0));
    w_regular = subs(subs(w{i}, dirac(1,t), sym(0)), dirac(t), sym(0));

    deltaH = simplify(h{i} - h_regular);
    deltaI = simplify(w{i} - w_regular);

    if ~isequal(deltaH, sym(0))
        noteH{i} = ['Добавка: ' char(simplify(deltaH/K*K_symbolic))];
    end
    if ~isequal(deltaI, sym(0))
        noteI{i} = ['Добавка: ' char(simplify(deltaI/K*K_symbolic))];
    end
end

colors = lines(numel(K_values));
for mode = 1:2
    figure('Color','w','Position',[50 50 1450 950]);
    tiledlayout(4,3,'TileSpacing','compact','Padding','compact');

    for i = 1:11
        nexttile; hold on; grid on;
        if mode == 1
            base = H{i}; notes = noteH;
        else
            base = I{i}; notes = noteI;
        end
        for j = 1:numel(K_values)
            plot(time, base*(K_values(j)/K), 'Color',colors(j,:), 'LineWidth',1.5);
        end
        title(sprintf('%d. %s',i,names{i}),'FontSize',9);
        xlabel('Время, с');
        if mode == 1, ylabel('h(t)'); else, ylabel('w(t)'); end
        if ~isempty(notes{i})
            text(0.03,0.90,notes{i},'Units','normalized','FontSize',7,'Color',[0.7 0.1 0.1],'Interpreter','none');
        end
        hold off;
    end

    nexttile; hold on; axis off;
    p = gobjects(1,numel(K_values));
    for j = 1:numel(K_values)
        p(j) = plot(nan,nan,'Color',colors(j,:),'LineWidth',2);
    end
    labels = arrayfun(@(k)sprintf('K = %.2f',k),K_values,'UniformOutput',false);
    legend(p,labels,'Location','northwest','FontSize',10);
    if mode == 1
        sgtitle('Переходные характеристики 11 звеньев: влияние K');
    else
        sgtitle('Импульсные характеристики 11 звеньев: влияние K');
    end
end

T_values = [T/5,T/2,T,2*T,5*T];
figure('Color','w'); hold on; grid on;
for Ti = T_values
    W4 = K/(Ti*s+1);
    h4 = ilaplace(W4*U_step, s, t);
    y = regular_values(h4, t, time);
    plot(time,y,'LineWidth',1.5,'DisplayName',sprintf('T = %.2f',Ti));
end
title('Звено №4: влияние постоянной времени T');
xlabel('Время, с'); ylabel('h(t)'); legend('show','Location','southeast');
hold off;

T1_values = [T1,T1/2,T1/6];
time6 = 0:0.01:40;
figure('Color','w'); hold on; grid on;
for T1i = T1_values
    W6 = K/(T2^2*s^2+T1i*s+1);
    h6 = ilaplace(W6*U_step, s, t);
    y = regular_values(h6, t, time6);
    plot(time6,y,'LineWidth',1.5,'DisplayName',sprintf('T1 = %.2f',T1i));
end
title('Звено №6: влияние T1 при постоянном T2');
xlabel('Время, с'); ylabel('h(t)'); legend('show','Location','northeast');
hold off;

G = {
    tf(K,[1 0]);
    tf([K 0],1);
    tf(K,1);
    tf(K,[T 1]);
    tf(K,[T2^2 T1 1]);
    tf(K,[T2^2 T1 1]);
    tf(K,[T 0 1]);
    tf(K,[T 1 0]);
    tf([K 0],[T 1]);
    tf([K*T K],1);
    tf([K*T K],[1 0])
};

fprintf('\n========== ЧИСЛЕННАЯ ПРОВЕРКА step / impulse ==========\n');
for i = 1:11
    fprintf('Звено %2d:\n', i);

    if contains(char(h{i}), 'dirac')
        fprintf('  h(t): пропуск (есть delta-составляющая)\n');
    else
        y = step(G{i}, time);
        err = max(abs(y(:).' - H{i}));
        if isfinite(err)
            fprintf('  h(t): max |h - step(G)| = %.3e\n', err);
        else
            fprintf('  h(t): результат не является конечным числом\n');
        end
    end

    if contains(char(w{i}), 'dirac')
        fprintf('  w(t): пропуск (есть delta-составляющая)\n');
    else
        y = impulse(G{i}, time);
        err = max(abs(y(:).' - I{i}));
        if isfinite(err)
            fprintf('  w(t): max |w - impulse(G)| = %.3e\n', err);
        else
            fprintf('  w(t): результат не является конечным числом\n');
        end
    end
end

function y = regular_values(expr, t, time)
    expr = subs(expr, dirac(1,t), sym(0));
    expr = subs(expr, dirac(t), sym(0));
    f = matlabFunction(expr, 'Vars', t);
    y = f(time);
    if isscalar(y)
        y = y * ones(size(time));
    end
    y = reshape(double(y), size(time));
end
