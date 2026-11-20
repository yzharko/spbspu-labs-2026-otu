

clc;
clear;
close all;

pkg load control;
pkg load symbolic;

syms s;

N = 100;
D = 5*s^4 + sym(1)/10*s^3 + 2*s^2 + 2*s + 1;

Wopen = tf(N, double(sym2poly(D)));
Wclosed = feedback(Wopen, 1);

Dopen = D;
Dclosed = D + N;

disp('Передаточная функция разомкнутой системы:');
disp(Wopen);
disp('Передаточная функция замкнутой системы:');
disp(Wclosed);
disp('Характеристический полином разомкнутой системы:');
disp(Dopen);
disp('Характеристический полином замкнутой системы:');
disp(Dclosed);

a = double(sym2poly(Dopen));
b = double(sym2poly(Dclosed));
disp('Коэффициэнты характеристического полинома разомкнутой системы:');
disp(a);
disp('Коэффициэнты характеристического полинома замкнутой системы:');
disp(b);

Hopen = [a(2), a(4), 0,    0;
         a(1), a(3), a(5), 0;
         0,    a(2), a(4), 0;
         0,    a(1), a(3), a(5)];

Hclosed = [b(2), b(4), 0,    0;
           b(1), b(3), b(5), 0;
           0,    b(2), b(4), 0;
           0,    b(1), b(3), b(5)];

disp('Матрица Гурвица разомкнутой системы:');
disp(Hopen);
disp('Матрица Гурвица замкнутой системы:');
disp(Hclosed);

disp('Определители разомкнутой системы:');
for i = 1:4
    disp(det(Hopen(1:i, 1:i)));
end

disp('Определители замкнутой системы:');
for i = 1:4
    disp(det(Hclosed(1:i, 1:i)));
end

Wopen_tf = tf(N, a);
Wclosed_tf = tf(N, b);
w = [0, logspace(-3, 2, 3000)];
[re, im] = nyquist(Wopen_tf, w);

figure(1);
plot(re(:), im(:), 'b', 'LineWidth', 1.5);
hold on;
plot(re(:), -im(:), 'r--', 'LineWidth', 1.5);
plot(-1, 0, 'kx', 'MarkerSize', 10, 'LineWidth', 2);
grid on;
axis equal;
xlabel('Re W(j\omega)');
ylabel('Im W(j\omega)');
title('Годограф Найквиста разомкнутой системы');
legend('\omega >= 0', '\omega <= 0', ...
       '(-1; 0)', 'Location', 'southoutside');
hold off;

figure(2);
plot(re(:), im(:), 'b', 'LineWidth', 1.5);
hold on;
plot(re(:), -im(:), 'r--', 'LineWidth', 1.5);
plot(-1, 0, 'kx', 'MarkerSize', 10, 'LineWidth', 2);
grid on;
axis equal;
xlim([-2 5]);
ylim([-3 3]);
xlabel('Re W(j\omega)');
ylabel('Im W(j\omega)');
title('Годограф Найквиста вблизи точки (-1; 0)');
legend('\omega >= 0', '\omega <= 0', ...
       '(-1; 0)', 'Location', 'southoutside');
hold off;

popen = pole(Wopen_tf);
P = sum(real(popen) > 0);
N_ccw = 0;
Z = P - N_ccw;
disp('Устойчивость замкнутой системы по Найквисту:');
disp(Z == 0);

w = linspace(0, 6, 3000);
D_mikhailov = polyval(b, 1i*w);
U = real(D_mikhailov);
V = imag(D_mikhailov);
figure(3);
plot(U, V, 'b', 'LineWidth', 1.5);
hold on;
plot(U(1), V(1), 'ko', 'MarkerFaceColor', 'k');
grid on;
xlabel('Re D(j\omega)');
ylabel('Im D(j\omega)');
title('Годограф Михайлова замкнутой системы');
legend('D(j\omega)', '\omega = 0', ...
       'Location', 'southoutside');
hold off;

figure(4);
plot(U, V, 'b', 'LineWidth', 1.5);
hold on;
plot(U(1), V(1), 'ko', 'MarkerFaceColor', 'k');
grid on;
xlim([99 180]);
ylim([-1 4]);
xlabel('Re D(j\omega)');
ylabel('Im D(j\omega)');
title('Начальный участок годографа Михайлова');
legend('D(j\omega)', '\omega = 0', ...
       'Location', 'southoutside');
hold off;

popen = pole(Wopen_tf);
pclosed = pole(Wclosed_tf);
disp('Полюса разомкнутой системы:');
disp(popen);
disp('Полюса замкнутой системы:');
disp(pclosed);

t = 0:0.01:12;
hopen = step(Wopen_tf, t);
hclosed = step(Wclosed_tf, t);
figure(5);
plot(t, hopen, 'b', 'LineWidth', 1.5);
grid on;
xlabel('Время, с');
ylabel('h(t)');
title('Переходная характеристика разомкнутой системы');
figure(6);
plot(t, hclosed, 'r', 'LineWidth', 1.5);
grid on;
xlabel('Время, с');
ylabel('h(t)');
title('Переходная характеристика замкнутой системы');

beta_open = abs(imag(popen(real(popen) > 0)));
beta_closed = abs(imag(pclosed(real(pclosed) > 0)));
Topen = 2*pi / beta_open(1);
Tclosed = 2*pi / beta_closed(1);
disp('Период колебаний разомкнутой системы:');
disp(Topen);
disp('Период колебаний замкнутой системы:');
disp(Tclosed);

