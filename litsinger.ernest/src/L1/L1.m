clc; clear; close all;

K  = 0.8;
T  = 0.5;
T1 = 0.7;
T2 = 0.5;
T1a = 1.5;

t_end = 10;
time = 0:0.01:t_end;

name = { ...
 'Интегрирующее', 'Дифференцирующее', 'Усилительное', ...
 'Апериодическое 1-го порядка', 'Апериодическое 2-го порядка', 'Колебательное', ...
 'Консервативное', 'Инерционное интегрирующее', 'Реальное дифференцирующее', ...
 'Форсирующее 1-го порядка', 'Изодромное'};
N = numel(name);

syms s t Ks real
Wsym = { Ks/s, Ks*s, Ks, Ks/(T*s+1), Ks/(T2^2*s^2+T1a*s+1), ...
         Ks/(T2^2*s^2+T1*s+1), Ks/(T*s^2+1), Ks/(s*(T*s+1)), ...
         Ks*s/(T*s+1), Ks*(T*s+1), Ks*(T*s+1)/s };

for i = 1:N
    h = ilaplace(Wsym{i}/s, s, t);
    w = ilaplace(Wsym{i},   s, t);
    fprintf('\n[%d] %s\n  h(t) = %s\n  w(t) = %s\n', ...
            i, name{i}, char(simplify(h)), char(simplify(w)));
end

hReg = logical([1 0 1 1 1 1 1 1 1 0 1]);
wReg = logical([1 0 0 1 1 1 1 1 0 0 0]);
hLbl = {'', 'h(t)=K\delta(t)', '', '', '', '', '', '', '', ...
        'h(t)=KT\delta(t)+K\cdot1(t)', ''};
wLbl = {'', 'w(t)=K\delta''(t)', 'w(t)=K\delta(t)', '', '', '', '', '', ...
        'w(t)=(K/T)\delta(t)-(K/T^2)e^{-t/T}', ...
        'w(t)=KT\delta''(t)+K\delta(t)', 'w(t)=KT\delta(t)+K\cdot1(t)'};

s_tf = tf('s');
W = { K/s_tf, K*s_tf, tf(K,1), K/(T*s_tf+1), K/(T2^2*s_tf^2+T1a*s_tf+1), ...
      K/(T2^2*s_tf^2+T1*s_tf+1), K/(T*s_tf^2+1), K/(s_tf*(T*s_tf+1)), ...
      K*s_tf/(T*s_tf+1), K*(T*s_tf+1), K*(T*s_tf+1)/s_tf };

figure('Name','Переходные характеристики h(t)','Position',[50 50 1200 720]);
for i = 1:N
    subplot(3,4,i); grid on;
    if hReg(i)
        [y,tt] = step(W{i}, time);
        plot(tt, y, 'LineWidth', 1.4);
    else
        axis([0 1 0 1]); text(0.5, 0.55, hLbl{i}, 'HorizontalAlignment','center');
        text(0.5, 0.35, '(график не строится)', 'HorizontalAlignment','center');
    end
    title(sprintf('%d. %s', i, name{i}), 'FontSize', 8);
    xlabel('t, с'); ylabel('h(t)');
end
sgtitle('Переходные характеристики типовых звеньев, K = 0.8');

figure('Name','Импульсные характеристики w(t)','Position',[50 50 1200 720]);
for i = 1:N
    subplot(3,4,i); grid on;
    if wReg(i)
        [y,tt] = impulse(W{i}, time);
        plot(tt, y, 'LineWidth', 1.4);
    else
        axis([0 1 0 1]); text(0.5, 0.55, wLbl{i}, 'HorizontalAlignment','center');
        text(0.5, 0.35, '(график не строится)', 'HorizontalAlignment','center');
    end
    title(sprintf('%d. %s', i, name{i}), 'FontSize', 8);
    xlabel('t, с'); ylabel('w(t)');
end
sgtitle('Импульсные характеристики типовых звеньев, K = 0.8');

figure('Name','Влияние K','Position',[80 80 700 500]); hold on; grid on;
for k = 1:10
    [y,tt] = step(k/(T2^2*s_tf^2+T1*s_tf+1), time);
    plot(tt, y, 'DisplayName', ['K = ', num2str(k)]);
end
title('Колебательное звено при K = 1..10');
xlabel('Время t, с'); ylabel('h(t)'); legend('show','Location','northwest'); hold off;

figure('Name','Влияние T','Position',[80 80 700 500]); hold on; grid on;
for Ti = [T/5, T/2, T, 2*T, 5*T]
    [y,tt] = step(K/(Ti*s_tf+1), time);
    plot(tt, y, 'DisplayName', ['T = ', num2str(Ti)]);
end
title('Апериодическое звено 1-го порядка при разных T');
xlabel('Время t, с'); ylabel('h(t)'); legend('show','Location','southeast'); hold off;

figure('Name','Влияние T1','Position',[80 80 700 500]); hold on; grid on;
for T1i = [0.7, 0.4, 0.2, 0.05]
    [y,tt] = step(K/(T2^2*s_tf^2+T1i*s_tf+1), time);
    plot(tt, y, 'DisplayName', ['T1 = ', num2str(T1i)]);
end
title('Колебательное звено при уменьшении T1 (T2 = 0.5)');
xlabel('Время t, с'); ylabel('h(t)'); legend('show','Location','northeast'); hold off;

h4 = K*(1 - exp(-time/T));
[y4,tt4] = step(W{4}, time);
figure('Name','Проверка h(t)','Position',[80 80 700 500]);
plot(time, h4, 'LineWidth', 2); hold on; grid on;
plot(tt4, y4, '--', 'LineWidth', 2);
title('Проверка переходной функции апериодического звена');
xlabel('Время t, с'); ylabel('h(t)');
legend('аналитически K(1-e^{-t/T})','step()','Location','southeast'); hold off;

