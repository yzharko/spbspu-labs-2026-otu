clc;
clear;
close all;

den = [0.05 0.1 1 1 1];
W = tf(1, den);
W_closed = feedback(W, 1);
denClosed = den + [0 0 0 0 1];

% Критерий Гурвица и корни характеристических полиномов
show_hurwitz('Разомкнутая система', den);
show_hurwitz('Замкнутая система', denClosed);

% Годограф Найквиста разомкнутой системы
omegaCross = sqrt(den(4)/den(2));
omega = unique([0 logspace(-3, 3, 2000) omegaCross]);
[re, im] = nyquist(W, omega);
re = squeeze(re);
im = squeeze(im);
response = re + 1i*im;
phase = unwrap(angle(1 + [flipud(conj(response)); response]));
turns = round(abs(phase(end) - phase(1))/(2*pi));
fprintf('\nНайквист: W(j*%.4f) = %.4f, число охватов (-1, 0) = %d\n', ...
    omegaCross, real(evalfr(W, 1i*omegaCross)), turns);

figure;
plot(re, im, 'b', re, -im, 'b');
hold on;
plot(-1, 0, 'r+', 'MarkerSize', 10);
grid on;
axis equal;
xlabel('Re W(j\omega)');
ylabel('Im W(j\omega)');
title('Годограф Найквиста разомкнутой системы');

% Годограф Михайлова замкнутой системы
omega = 0:0.01:6;
Djw = polyval(denClosed, 1i*omega);
highFreq = logspace(log10(6), 3, 1000);
phase = unwrap(angle([Djw, polyval(denClosed, 1i*highFreq)]));
fprintf('Михайлов: при 1000 рад/с поворот аргумента %.2f град\n', ...
    rad2deg(phase(end) - phase(1)));

figure;
plot(real(Djw), imag(Djw));
hold on;
xline(0, 'k:');
yline(0, 'k:');
grid on;
xlabel('Re D(j\omega)');
ylabel('Im D(j\omega)');
title('Годограф Михайлова замкнутой системы');

% Переходные характеристики
t = (0:0.02:20).';
yOpen = step(W, t);
yClosed = step(W_closed, t);
fprintf('Установившиеся значения: %.2f и %.2f\n', dcgain(W), dcgain(W_closed));

figure;
plot(t, yOpen, t, yClosed);
grid on;
xlabel('t, с');
ylabel('h(t)');
title('Отклики на единичный скачок');
legend('Разомкнутая система', 'Замкнутая система', 'Location', 'southeast');

function show_hurwitz(name, a)
    H = [a(2) a(4) 0    0;
         a(1) a(3) a(5) 0;
         0    a(2) a(4) 0;
         0    a(1) a(3) a(5)];
    delta = zeros(1, 4);
    for k = 1:4
        delta(k) = det(H(1:k, 1:k));
    end
    fprintf('\n%s\n', name);
    fprintf('Определители Гурвица: %.4f %.4f %.4f %.4f\n', delta);
    disp('Корни характеристического полинома:');
    disp(roots(a));
end
