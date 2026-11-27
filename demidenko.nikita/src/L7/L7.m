A = [1.0 -15; 0.5 -2.0];
B = [2.0; 0.0];
C = [0.5 1.0];

M = [2.0 3.0; 0.0 5.0];
A_new = M \ A * M;
B_new = M \ B;
C_new = C * M;
show_system("Новый базис", A_new, B_new, C_new);
% A: [0.7 -33.45;0.2 -1.7]
% B: [1;0]
% C: [1 6.5]

syms s
W = C * ((eye(2) * s - A) \ B);
W = simplifyFraction(W, "Expand", true);
disp(W); % (2*s + 6)/(2*s^2 + 2*s + 11)
[num, den] = numden(W);
num = sym2poly(num); den = sym2poly(den);
[auto_num, auto_den] = ss2tf(A, B, C, 0);

fprintf("num: %s, den: %s\n", mat2str(num), mat2str(den));           
fprintf("num: %s, den: %s\n", mat2str(auto_num), mat2str(auto_den));
% num:   [2 6], den:  [2 2 11]
% num: [0 1 3], den: [1 1 5.5]

init_system = ss(A, B, C, 0);
new_system = ss(A_new, B_new, C_new, 0);
tf_system = tf(num, den);

figure("Name", "Сравнение");
hold on;
[y,t] = step(init_system); plot(t, y, "LineWidth", 6);
[y,t] = step(new_system);  plot(t, y, "LineWidth", 4);
[y,t] = step(tf_system);   plot(t, y, "LineWidth", 2);
hold off;
legend({
    "Изначальная"
    "Новый базис"
    "Передаточная функция"
});

ctrb_mat = ctrb(A, B);
fprintf("Матрица управляемости:\n");
disp(ctrb_mat);
% 2     2
% 0     1
fprintf("Ранг: %d\n", rank(ctrb_mat)); % 2 - управляема

obsv_mat = obsv(A, C);
fprintf("Матрица наблюдаемости:\n");
disp(obsv_mat);
% 0.5000    1.0000
% 1.0000   -9.5000
fprintf("Ранг: %d\n", rank(obsv_mat)); % 2 - наблюдаема

num = num / den(1);
den = den / den(1);

A_ctrb = [0 1; -den(3), -den(2)];
B_ctrb = [0; 1];
C_ctrb = flip(num);
ctrb_system = ss(A_ctrb, B_ctrb, C_ctrb, 0);
show_system("Управляемая", A_ctrb, B_ctrb, C_ctrb);
% A: [0 1;-5.5 -1]
% B: [0;1]
% C: [3 1]

A_obsv = A_ctrb';
B_obsv = C_ctrb';
C_obsv = B_ctrb';
obsv_system = ss(A_obsv, B_obsv, C_obsv, 0);
show_system("Наблюдаемая", A_obsv, B_obsv, C_obsv);
% A: [0 -5.5;1 -1]
% B: [3;1]
% C: [0 1]

figure("Name", "Сравнение2");
hold on;
[y,t] = step(init_system); plot(t, y, "LineWidth", 6);
[y,t] = step(ctrb_system); plot(t, y, "LineWidth", 4);
[y,t] = step(obsv_system); plot(t, y, "LineWidth", 2);
hold off;
legend({
    "Изначальная"
    "Управляемая"
    "Наблюдаемая"
});

function show_system(name, A, B, C)
    fprintf("%s:\n A: %s\n B: %s\n C: %s\n\n", ...
        name, mat2str(A), mat2str(B), mat2str(C));
end
