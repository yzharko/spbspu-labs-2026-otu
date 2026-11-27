clear;
clc;
close all;

n = 2;
A  = [1, -15;
      1, -3];
B  = [2; 0];
CT = [2; 1];
C  = CT';      
M  = [5, 0;
      6, 2];

M_inv = inv(M);
A_new = M_inv * A * M;
B_new = M_inv * B;
C_new = C * M;

[num, den] = ss2tf(A, B, C, 0);
W_tf = tf(num, den);

[ya_num, ya_dem] = ss2tf(A_new, B_new, C_new, 0);
ya_W_tf = tf(ya_num, ya_dem);

t = 0:0.01:10;
sys_old = ss(A, B, C, 0);
sys_new = ss(A_new, B_new, C_new, 0);

[y_old, ~] = step(sys_old, t);
[y_new, ~] = step(sys_new, t);
[y_tf, ~] = step(W_tf, t);

fig1 = figure('Name','Переходные процессы (3 способа)','Color','w');
plot(t, y_old, 'r-', 'LineWidth', 1.5); hold on;
plot(t, y_new, 'b--', 'LineWidth', 1.5);
plot(t, y_tf, 'g:', 'LineWidth', 2);
grid on; xlabel('t, с'); ylabel('y(t)');
title('Переходной процесс при единичном ступенчатом воздействии');
legend('По исходной модели','По модели в новом базисе','По передаточной функции W(s)','Location','southeast');
saveas(fig1, 'Fig1_3methods.png');

b1 = num(2); b0 = num(3);
a1 = den(2); a0 = den(3);

Ac = [0 -a0; 1 -a1];
Bc = [1; 0];
Cc = [b1, b0 - a1*b1];

Ao = [0 1; -a0 -a1];
Bo = [b1; b0 - a1*b1];
Co = [1 0];

sys_c = ss(Ac, Bc, Cc, 0);
sys_o = ss(Ao, Bo, Co, 0);
y_c = step(sys_c, t);
y_o = step(sys_o, t);

fig2 = figure('Name','Канонические формы','Color','w');
plot(t, y_old, 'r-', 'LineWidth', 1.5); hold on;
plot(t, y_c, 'b--', 'LineWidth', 1.5);
plot(t, y_o, 'g:', 'LineWidth', 2);
grid on; xlabel('t, с'); ylabel('y(t)');
title('Сравнение переходных характеристик');
legend('Исходная','Канонич. управляемая','Канонич. наблюдаемая','Location','southeast');
saveas(fig2, 'Fig2_canonical.png');

