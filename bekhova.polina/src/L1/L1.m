clc; clear; close all;
format compact;

K = 3.8; T = 1.1; T1 = 1.5; T2 = 1.6;
syms s t

T1_ap = 3.5;

U_step = simplify(laplace(heaviside(t), t, s));    
U_impulse = simplify(laplace(dirac(t), t, s));

fprintf('U_step(s) = '); disp(U_step);
fprintf('U_impulse(s) = '); disp(U_impulse);

names = {'Интегрирующее', 'Дифференцирующее', 'Усилительное', ...
         'Апериодическое 1-го порядка', 'Апериодическое 2-го порядка', ...
         'Колебательное', 'Консервативное', 'Инерционное интегрирующее', ...
         'Реальное дифференцирующее', 'Форсирующее 1-го порядка', 'Изодромное'};

N = numel(names);

W = {K/s;                   
     K*s;                  
     K;    
     K/(T*s+1);       
     K/(T2^2*s^2+T1_ap*s+1); 
     K/(T2^2*s^2+T1*s+1); 
     K/(T*s^2+1);      
     K/(s*(T*s+1));   
     K*s/(T*s+1);  
     K*(T*s+1);   
     K*(T*s+1)/s};   

h = cell(N,1); w = cell(N,1);
for i = 1:N
    Y_step = simplify(W{i} * U_step);
    h{i} = simplify(ilaplace(Y_step, s, t));
    Y_impulse = simplify(W{i} * U_impulse);
    w{i} = simplify(ilaplace(Y_impulse, s, t));

    fprintf('\nЗвено %d: %s\n', i, names{i});
    fprintf('W(s) = '); disp(W{i});
    fprintf('h(t) = '); disp(h{i});
    fprintf('w(t) = '); disp(w{i});
end

fprintf('\nАналитическая проверка\n');
for i = 1:N
    Y_step = simplify(W{i} * U_step);
    Y_impulse = simplify(W{i} * U_impulse);
    check_h = simplify(laplace(h{i}, t, s) - Y_step);
    check_w = simplify(laplace(w{i}, t, s) - Y_impulse);
    fprintf('Звено %d: h(t) -> Y_step = %s; w(t) -> Y_impulse = %s\n', ...
        i, char(check_h), char(check_w));
end

time = 0:0.01:12;
Hplot = zeros(N, numel(time));
Wplot = zeros(N, numel(time));

Hplot(1,:) = K*time;
Hplot(3,:) = K*ones(size(time));
Hplot(4,:) = K*(1-exp(-time/T));

r1 = (-T1_ap+sqrt(T1_ap^2-4*T2^2))/(2*T2^2);
r2 = (-T1_ap-sqrt(T1_ap^2-4*T2^2))/(2*T2^2);
h5 = K*(1+(r2*exp(r1*time)-r1*exp(r2*time))/(r1-r2));
w5 = K/T2^2*(exp(r1*time)-exp(r2*time))/(r1-r2);
Hplot(5,:) = real(h5); Wplot(5,:) = real(w5);

alpha = T1/(2*T2^2); 
beta = sqrt(4*T2^2-T1^2)/(2*T2^2);
h6 = K*(1-exp(-alpha*time).*(cos(beta*time)+(alpha/beta)*sin(beta*time)));
w6 = K/(T2^2*beta)*exp(-alpha*time).*sin(beta*time);
Hplot(6,:) = real(h6); Wplot(6,:) = real(w6);

Hplot(7,:) = K*(1-cos(time/sqrt(T)));
Hplot(8,:) = K*(time-T+T*exp(-time/T));
Hplot(9,:) = (K/T)*exp(-time/T);
Hplot(10,:) = K*ones(size(time));
Hplot(N,:) = K*(time+T);

Wplot(1,:) = K*ones(size(time));
Wplot(4,:) = (K/T)*exp(-time/T);
Wplot(7,:) = (K/sqrt(T))*sin(time/sqrt(T));
Wplot(8,:) = K*(1-exp(-time/T));
Wplot(9,:) = -(K/T^2)*exp(-time/T);
Wplot(N,:) = K*ones(size(time));

figure('Name','Переходные характеристики h(t)','Color','w');
cols = ceil(sqrt(N)); rows = ceil(N/cols);
tiledlayout(rows, cols,'TileSpacing','compact');
for i = 1:N
    nexttile; plot(time,Hplot(i,:),'LineWidth',1.4); grid on;
    xlabel('t, с'); ylabel('h(t)');
    title(sprintf('%d. %s',i,names{i}),'FontSize',8);
end

figure('Name','Импульсные характеристики w(t)','Color','w');
cols = ceil(sqrt(N)); rows = ceil(N/cols);
tiledlayout(rows, cols,'TileSpacing','compact');
for i = 1:N
    nexttile; plot(time,Wplot(i,:),'LineWidth',1.4); grid on;
    xlabel('t, с'); ylabel('w(t)');
    title(sprintf('%d. %s',i,names{i}),'FontSize',8);
end

s_tf = tf('s');
W_tf = {K/s_tf; K*s_tf; K; K/(T*s_tf+1); K/(T2^2*s_tf^2+T1_ap*s_tf+1);
        K/(T2^2*s_tf^2+T1*s_tf+1); K/(T*s_tf^2+1); K/(s_tf*(T*s_tf+1));
        K*s_tf/(T*s_tf+1); K*(T*s_tf+1); K*(T*s_tf+1)/s_tf};

figure('Name','Проверка CST: step','Color','w'); hold on; grid on;
for i = [1 3 4 5 6 7 8 9 N]
    step(W_tf{i}, time);
end
title('Проверка переходных характеристик (CST)');
xlabel('t, с'); ylabel('h(t)');
legend(string([1 3 4 5 6 7 8 9 N])); hold off;

figure('Name','Проверка CST: impulse','Color','w'); hold on; grid on;
for i = [1 4 5 6 7 8]
    impulse(W_tf{i}, time);
end
title('Проверка импульсных характеристик (CST)');
xlabel('t, с'); ylabel('w(t)');
legend(string([1 4 5 6 7 8])); hold off;

T_vals = [T/5, T/2, T, 2*T, 5*T];
figure('Name','Звено №4. влияние T','Color','w'); hold on; grid on;
for Tv = T_vals
    h4 = K*(1-exp(-time/Tv));
    plot(time, h4, 'LineWidth', 1.5, 'DisplayName', sprintf('T=%.3f', Tv));
end
title('Апериодическое звено №4: влияние T');
xlabel('t, с'); ylabel('h(t)'); legend('Location','southeast'); hold off;

fprintf('\nОценка скорости протекания процесса:\n');
for Tv = T_vals
    fprintf('T = %.2f с: t95 ≈ %.2f с\n', Tv, 3*Tv);
end

T1_vals = [1.5 1.0 0.5 0.2 0.1 0];
figure('Name','Звено №6. влияние T1','Color','w'); hold on; grid on;
for T1v = T1_vals
    if T1v == 0
        h6 = K*(1-cos(time/T2));
    else
        a = T1v/(2*T2^2); b = sqrt(4*T2^2-T1v^2)/(2*T2^2);
        h6 = K*(1-exp(-a*time).*(cos(b*time)+(a/b)*sin(b*time)));
    end
    plot(time, h6, 'LineWidth', 1.3, 'DisplayName', sprintf('T1=%.1f', T1v));
end
title('Колебательное звено №6: влияние T1 (T2=const)');
xlabel('t, с'); ylabel('h(t)'); legend('Location','eastoutside'); hold off;

K_vals = [0.5 1 1.5 3 5];
figure('Name','Влияние K','Color','w');
tiledlayout(2,2,'TileSpacing','compact');

nexttile; hold on; grid on;
for kv = K_vals, plot(time, kv*time, 'LineWidth',1.3, 'DisplayName', sprintf('K=%.1f',kv)); end
title('№1 Интегрирующее'); xlabel('t, с'); ylabel('h(t)');
legend('Location','northwest'); hold off;

nexttile; hold on; grid on;
for kv = K_vals, plot(time, kv*ones(size(time)), 'LineWidth',1.3, 'DisplayName', sprintf('K=%.1f',kv)); end
title('№3 Усилительное'); xlabel('t, с'); ylabel('h(t)');
legend('Location','best'); hold off;

nexttile; hold on; grid on;
for kv = K_vals, plot(time, kv*(1-exp(-time/T)), 'LineWidth',1.3, 'DisplayName', sprintf('K=%.1f',kv)); end
title('№4 Апериодическое'); xlabel('t, с'); ylabel('h(t)');
legend('Location','southeast'); hold off;

nexttile; hold on; grid on;
a = T1/(2*T2^2); b = sqrt(4*T2^2-T1^2)/(2*T2^2);
for kv = K_vals
    h6 = kv*(1-exp(-a*time).*(cos(b*time)+(a/b)*sin(b*time)));
    plot(time, h6, 'LineWidth',1.3, 'DisplayName', sprintf('K=%.1f',kv));
end
title('№6 Колебательное'); xlabel('t, с'); ylabel('h(t)');
legend('Location','southeast'); hold off;

sgtitle('Влияние коэффициента K на переходные характеристики');
