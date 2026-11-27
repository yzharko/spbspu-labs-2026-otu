clear; clc; close all;

s = tf("s");
Go = 2 / (s * (s + 1) * (s + 8));
disp_margin(Go); % Gm = 36, Pm = 74.6126

Go_closed = feedback(Go, 1);
disp(mat2str(gurv(Go_closed.Denominator))); % [9 70 140] Устойчива
disp(stepinfo(Go_closed));
%      RiseTime: 6.4125
% TransientTime: 11.0598
%  SettlingTime: 11.0598
%   SettlingMin: 0.9020
%   SettlingMax: 0.9996
%     Overshoot: 0
%    Undershoot: 0
%          Peak: 0.9996
%      PeakTime: 18.0460

t_range = 0:0.1:50;
K_range = [1, 20, 35, 36, 37]; 
kvariation(Go, t_range, K_range);
K_crit = 36; % незатухающие колебания

figure("Name", "Определение периода колебаний");
G_osc = feedback(K_crit * Go, 1);
t_range = 0:0.01:10;
y_range = step(G_osc, t_range);
plot(t_range, y_range); hold on;
plot([0, 10], [1, 1], 'k');
d = y_range - 1;
intersect = find(d(1:end-1) .* d(2:end) <= 0)';
for i = intersect
    t = t_range(i);
    plot([t, t], [0, 1], 'k--');
end
xticks(0:0.2:10);
hold off;
T_crit = 2.2; % По графику

K_r = 0.6 * K_crit;
T_i = 0.5 * T_crit;
T_d = 0.125 * T_crit;

W_reg = K_r * (1 + 1 / (T_i * s) + T_d * s);
Go_reg = W_reg * Go;
disp_margin(Go_reg); %Gm = 0.0256799, Pm = 19.8915
Go_reg_closed = feedback(Go_reg, 1);
figure("Name", "Отрегулированный");
step(Go_reg_closed);
disp(stepinfo(Go_reg_closed));
%      RiseTime: 0.4492
% TransientTime: 9.3717
%  SettlingTime: 9.3717
%   SettlingMin: 0.6411
%   SettlingMax: 1.6394
%     Overshoot: 63.9361
%    Undershoot: 0
%          Peak: 1.6394
%      PeakTime: 1.3065

C = pidtune(Go, "PID");
disp(C); % Kp = 8.9161; Ki = 1.5697; Kd = 7.1841
Go_auto = C * Go;
disp_margin(Go_auto); % Gm = 0, Pm = 69.9597
Go_auto_closed = feedback(Go_auto, 1);
disp(stepinfo(Go_auto_closed));
%      RiseTime: 0.7648
% TransientTime: 8.7567
%  SettlingTime: 8.7567
%   SettlingMin: 0.9074
%   SettlingMax: 1.0945
%     Overshoot: 9.4505
%    Undershoot: 0
%          Peak: 1.0945
%      PeakTime: 2.1906

figure("Name", "Зиглер-Николс и автоматический");
step(Go_reg_closed); hold on;
step(Go_auto_closed); hold off;
legend("З-Н", "Авто");

function disp_margin(sys)
    [Gm, Pm] = margin(sys);
    fprintf("Gm = %g, Pm = %g\n", Gm, Pm);
end
function kvariation(sys, t_range, K_range)
    figure("Name", "Варьирование К");
    hold on;
    ylim([-1, 3]);
    for K = K_range
        closed = feedback(K * sys, 1);
        [y_range, ~] = step(closed, t_range);
        plot(t_range, y_range)
    end
    labels = arrayfun(@(K) strcat("K = ", num2str(K)), K_range);
    legend(labels);
    hold off;
end
function g = gurv(d)
    a = cell2mat(d);
    gurvMatrix = [
        a(2), a(4),    0
        a(1), a(3),    0
           0, a(2), a(4)
    ];
    first = a(2);
    second = det(gurvMatrix(1:2, 1:2));
    third = second * a(4);
    g = [first, second, third];
end
