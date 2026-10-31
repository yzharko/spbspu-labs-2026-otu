clear; clc; close all;
K = 2.6;
T = 0.3;
T1 = 0.6;
T2 = 0.5;
T1a = 1.125;
dt = 0.005;
time = 0:dt:10;
syms s t

names = { ...
    'Интегрирующее', 'Идеальное дифференцирующее', 'Усилительное', ...
    'Апериодическое 1-го порядка', 'Апериодическое 2-го порядка', ...
    'Колебательное', 'Консервативное', 'Инерционное интегрирующее', ...
    'Реальное дифференцирующее', 'Форсирующее 1-го порядка', 'Изодромное'};

W = {K/s, K*s, sym(K), K/(T*s+1), ...
    K/(T2^2*s^2+T1a*s+1), K/(T2^2*s^2+T1*s+1), ...
    K/(T*s^2+1), K/(s*(T*s+1)), K*s/(T*s+1), ...
    K*(T*s+1), K*(T*s+1)/s};

%% 1–6. Прямое и обратное преобразования Лапласа
U_step = laplace(heaviside(t), t, s);   % 1/s
U_impulse = laplace(dirac(t), t, s);  % 1
fprintf('U_step(s) = %s; U_impulse(s) = %s\n', char(U_step), char(U_impulse));
h = cell(1,11); w = cell(1,11);
for i = 1:11
    Y_step = simplify(W{i}*U_step);
    Y_impulse = simplify(W{i}*U_impulse);
    h{i} = simplify(ilaplace(Y_step, s, t));
    w{i} = simplify(ilaplace(Y_impulse, s, t));
    fprintf('\n%2d. %s\nW(s)=%s\nY_step(s)=%s\nh(t)=%s\n', ...
        i, names{i}, char(W{i}), char(Y_step), char(h{i}));
    fprintf('Y_impulse(s)=%s\nw(t)=%s\n', char(Y_impulse), char(w{i}));
end

%% 7–8. Переходные и импульсные характеристики
for input_kind = 1:2
    figure('Name', sprintf('ЛР1: вход %d', input_kind), 'Color', 'w');
    tiledlayout(4,3,'TileSpacing','compact','Padding','compact');
    for i = 1:11
        nexttile;
        y = regular_response(i,time,K,T,T1,T2,T1a,input_kind);
        plot(time,y,'b','LineWidth',1.3); grid on;
        xlim([0 10]);
        title(sprintf('%d. %s',i,names{i}),'FontSize',9);
        xlabel('t, c');
        if input_kind == 1, ylabel('h(t)'); else, ylabel('w(t)'); end
        singular = singular_label(i,input_kind,K,T);
        if ~isempty(singular)
            yl = ylim; dy = max(yl(2)-yl(1),1);
            text(0.2,yl(2)-0.12*dy,singular,'Color',[0.7 0 0], ...
                'FontSize',8,'Interpreter','tex');
        end
    end
    drawnow;
end

%% Влияние коэффициента K
s_tf = tf('s');
K_values = [K/2, K, 1.5*K, 2*K];
figure('Name','Влияние K на колебательное звено','Color','w');
hold on; grid on;
for k = K_values
    [y_k,t_k] = step(k/(T2^2*s_tf^2+T1*s_tf+1),time);
    plot(t_k(:),y_k(:),'LineWidth',1.3, ...
        'DisplayName',sprintf('K = %.2f',k));
end
title('Колебательное звено: изменение K');
xlabel('t, c'); ylabel('h(t)'); xlim([0 10]);
legend('show','Location','best');
drawnow;

%% Апериодическое звено: T/5, T/2, T, 2T, 5T
T_factors = [0.2, 0.5, 1, 2, 5];
figure('Name','Влияние T на апериодическое звено','Color','w');
hold on; grid on;
for factor = T_factors
    T_current = T*factor;
    [y_T,t_T] = step(K/(T_current*s_tf+1),time);
    plot(t_T(:),y_T(:),'LineWidth',1.3, ...
        'DisplayName',sprintf('T = %.2f c',T_current));
end
title('Апериодическое звено 1-го порядка: изменение T');
xlabel('t, c'); ylabel('h(t)'); xlim([0 10]);
legend('show','Location','best');
drawnow;

%% Колебательное звено: фиксируем T2, уменьшаем T1
T1_values = [T1, 0.45, 0.30, 0.20, 0.10];
time_T1 = 0:dt:30;
time_metrics = 0:dt:60;
metrics = zeros(numel(T1_values),4); 
figure('Name','Влияние T1 на колебательное звено','Color','w');
hold on; grid on;
for j = 1:numel(T1_values)
    T1_current = T1_values(j);
    system = K/(T2^2*s_tf^2+T1_current*s_tf+1);
    [y_plot,t_plot] = step(system,time_T1);
    plot(t_plot(:),y_plot(:),'LineWidth',1.3, ...
        'DisplayName',sprintf('T1 = %.2f c',T1_current));

    [y_metric,t_metric] = step(system,time_metrics);
    y_metric = y_metric(:);
    t_metric = t_metric(:);
    overshoot = max(0,(max(y_metric)-K)/K*100);
    last_outside = find(abs(y_metric-K)>0.02*K,1,'last');
    if isempty(last_outside)
        settling_time = 0;
    elseif last_outside < numel(t_metric)
        settling_time = t_metric(last_outside+1);
    else
        settling_time = NaN;
    end
    metrics(j,:) = [T1_current, T1_current/(2*T2), overshoot, settling_time];
end
title('Колебательное звено: уменьшение T1 при T2 = const');
xlabel('t, c'); ylabel('h(t)'); xlim([0 30]);
legend('show','Location','best');
drawnow;
fprintf('\nЗадание со звёздочкой: T2 = %.2f c, полоса установления 2%%\n',T2);
fprintf('    T1       zeta     overshoot, %%    settling, c\n');
for j = 1:size(metrics,1)
    fprintf('%7.3f    %7.3f     %10.2f      %10.3f\n',metrics(j,:));
end

function y = regular_response(i,t,K,T,T1,T2,T1a,input_kind)
    a = T1/(2*T2^2);
    b = sqrt(1-(T1/(2*T2))^2)/T2;
    discriminant = sqrt(T1a^2-4*T2^2);
    lambda1 = (T1a-discriminant)/(2*T2^2);
    lambda2 = (T1a+discriminant)/(2*T2^2);
    zero = zeros(size(t));
    one = ones(size(t));
    e = exp(-t/T);
    e2 = exp(-a*t);
    aperiodic_step = K*(1-(lambda2*exp(-lambda1*t) ...
        -lambda1*exp(-lambda2*t))/(lambda2-lambda1));
    aperiodic_impulse = K/(T2^2*(lambda2-lambda1)) ...
        *(exp(-lambda1*t)-exp(-lambda2*t));
    oscillatory_step = K*(1-e2.*(cos(b*t)+(a/b)*sin(b*t)));
    oscillatory_impulse = K/(T2^2*b)*e2.*sin(b*t);

    step_curves = {K*t, zero, K*one, K*(1-e), ...
        aperiodic_step, oscillatory_step, K*(1-cos(t/sqrt(T))), ...
        K*(t-T*(1-e)), K/T*e, K*one, K*(t+T)};
    impulse_curves = {K*one, zero, zero, K/T*e, ...
        aperiodic_impulse, oscillatory_impulse, ...
        K/sqrt(T)*sin(t/sqrt(T)), K*(1-e), ...
        -K/T^2*e, zero, K*one};

    if input_kind == 1
        y = step_curves{i};
    else
        y = impulse_curves{i};
    end
end

function label = singular_label(i,input_kind,K,T)
    step_labels = repmat({''},1,11);
    step_labels{2} = sprintf('В t=0: %.2g delta(t)',K);
    step_labels{10} = sprintf('В t=0: %.2g delta(t)',K*T);

    impulse_labels = repmat({''},1,11);
    impulse_labels{2} = sprintf('В t=0: %.2g delta''(t)',K);
    impulse_labels{3} = sprintf('В t=0: %.2g delta(t)',K);
    impulse_labels{9} = sprintf('В t=0: %.3g delta(t)',K/T);
    impulse_labels{10} = sprintf('В t=0: %.2g delta(t) + %.2g delta''(t)',K,K*T);
    impulse_labels{11} = sprintf('В t=0: %.2g delta(t)',K*T);

    if input_kind == 1
        label = step_labels{i};
    else
        label = impulse_labels{i};
    end
end
