clear; clc; close all;

K = 2.1;
T = 0.9;
T1 = 1.0;
T2 = 1.1;

T1_ap = 3.0; %для звена 5

t = 0:0.01:20;
z = zeros(size(t));
E = exp(-t/T);

D = sqrt(T1_ap^2 - 4*T2^2);
ta = (T1_ap + D)/2;
tb = (T1_ap - D)/2;
Ea = exp(-t/ta);
Eb = exp(-t/tb);

a = T1/(2*T2^2);
wd = sqrt(1/T2^2 - a^2);
Ed = exp(-a*t);
w0 = 1/sqrt(T);

links = {
    struct('name', 'Интегрирующее', 'h', K*t, 'w', K+z, 'delta', [0 0], 'dDelta', [0 0])
    struct('name', 'Идеальное дифференцирующее', 'h', z, 'w', z, 'delta', [K 0], 'dDelta', [0 K])
    struct('name', 'Усилительное', 'h', K+z, 'w', z, 'delta', [0 K], 'dDelta', [0 0])
    struct('name', 'Апериодическое 1-го порядка', 'h', K*(1-E), 'w', K/T*E, 'delta', [0 0], 'dDelta', [0 0])
    struct('name', 'Апериодическое 2-го порядка', 'h', K*(1-(ta*Ea-tb*Eb)/(ta-tb)), 'w', K*(Ea-Eb)/(ta-tb), 'delta', [0 0], 'dDelta', [0 0])
    struct('name', 'Колебательное', 'h', K*(1-Ed.*(cos(wd*t)+a/wd*sin(wd*t))), 'w', K/(T2^2*wd)*Ed.*sin(wd*t), 'delta', [0 0], 'dDelta', [0 0])
    struct('name', 'Консервативное', 'h', K*(1-cos(w0*t)), 'w', K*w0*sin(w0*t), 'delta', [0 0], 'dDelta', [0 0])
    struct('name', 'Инерционное интегрирующее', 'h', K*(t-T*(1-E)), 'w', K*(1-E), 'delta', [0 0], 'dDelta', [0 0])
    struct('name', 'Реальное дифференцирующее', 'h', K/T*E, 'w', -K/T^2*E, 'delta', [0 K/T], 'dDelta', [0 0])
    struct('name', 'Форсирующее', 'h', K+z, 'w', z, 'delta', [K*T K], 'dDelta', [0 K*T])
    struct('name', 'Изодромное', 'h', K*T+K*t, 'w', K+z, 'delta', [0 K*T], 'dDelta', [0 0])
};

fields = {'h', 'w'};
labels = {'h(t)', 'w(t)'};
figureNames = {'Реакция на единичный скачок', 'Реакция на единичный импульс'};

for j = 1:2
    figure('Name', figureNames{j}, 'NumberTitle', 'off');
    tiledlayout(4, 3);
    for n = 1:numel(links)
        link = links{n};
        nexttile;
        plot([-1 0 t], [0 0 link.(fields{j})], 'LineWidth', 1.2);
        grid on;
        xlim([-1 t(end)]);
        xlabel('t, с');
        ylabel(labels{j});
        title(sprintf('%d. %s', n, link.name));

        terms = {};
        if link.delta(j) ~= 0
            terms{end+1} = sprintf('%g\\delta(t)', link.delta(j));
        end
        if link.dDelta(j) ~= 0
            terms{end+1} = sprintf('%g\\delta''(t)', link.dDelta(j));
        end
        if ~isempty(terms)
            subtitle(['В нуле: ' strjoin(terms, ' + ')]);
        end
    end
    sgtitle(figureNames{j});
end


figure('Name', 'Влияние K', 'NumberTitle', 'off');
tiledlayout(1, 2);
scale = [0.5 1 2];
for j = 1:2
    nexttile;
    y = links{6}.(fields{j});
    plot(t, y(:) * scale, 'LineWidth', 1.2);
    grid on;
    xlabel('t, с');
    ylabel(labels{j});
    title('№6. Изменение K');
    legend(compose('K = %g', K*scale), 'Location', 'best');
end

T_test = T*[1/5 1/2 1 2 5];
t4 = 0:0.01:5*max(T_test);
E4 = exp(-t4(:)./T_test);
h4 = K*(1 - E4);
w4 = E4.*(K./T_test);

figure('Name', '№4. Влияние T', 'NumberTitle', 'off');
tiledlayout(1, 2);
curves4 = {h4, w4};
for j = 1:2
    nexttile;
    plot(t4, curves4{j}, 'LineWidth', 1.2);
    grid on;
    xlabel('t, с');
    ylabel(labels{j});
    legend(compose('T = %g с', T_test), 'Location', 'best');
    if j == 1
        yline(0.95*K, ':', '95% K', 'HandleVisibility', 'off');
    end
end
sgtitle('№4. Влияние постоянной времени T');

fprintf('\n№4. Апериодическое звено: время установления по зоне 5%%\n');
fprintf('%8s %12s %12s\n', 'T, с', 'tуст, с', 'w(0+)');
for Ti = T_test
    fprintf('%8.2f %12.3f %12.3f\n', Ti, Ti*log(20), K/Ti);
end

T1_test = T1*[1 1/2 1/4 1/10];
a6 = T1_test/(2*T2^2);
wd6 = sqrt(1/T2^2 - a6.^2);
t6 = linspace(0, 8/min(a6), 30001).';
E6 = exp(-t6.*a6);
S6 = sin(t6.*wd6);
error6 = E6.*(cos(t6.*wd6) + S6.*(a6./wd6));
h6 = K*(1 - error6);
w6 = E6.*S6.*(K./(T2^2*wd6));

figure('Name', '№6. Влияние T1', 'NumberTitle', 'off');
tiledlayout(1, 2);
curves6 = {h6, w6};
for j = 1:2
    nexttile;
    plot(t6, curves6{j}, 'LineWidth', 1.2);
    grid on;
    xlabel('t, с');
    ylabel(labels{j});
    legend(compose('T_1 = %g с', T1_test), 'Location', 'best');
    if j == 1
        yline(0.95*K, ':', 'HandleVisibility', 'off');
        yline(1.05*K, ':', 'HandleVisibility', 'off');
    end
end
sgtitle('№6. Ослабление затухания при фиксированном T_2');

fprintf('\n№6. Колебательное звено: перерегулирование и зона 5%%\n');
fprintf('%8s %10s %16s %12s\n', 'T1, с', 'zeta', 'Перерегул., %', 'tуст, с');
for n = 1:numel(T1_test)
    error_abs = abs(error6(:,n));
    last = find(error_abs > 0.05, 1, 'last');
    part = (error_abs(last) - 0.05)/(error_abs(last) - error_abs(last+1));
    settling = t6(last) + part*(t6(last+1) - t6(last));
    overshoot = 100*exp(-pi*a6(n)/wd6(n));
    fprintf('%8.2f %10.3f %16.2f %12.3f\n', ...
        T1_test(n), T1_test(n)/(2*T2), overshoot, settling);
end

