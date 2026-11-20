clc;
clear;
close all;

s = tf('s');
Wopen = 1/(0.1*s^3 + 0.1*s^2 + s + 1);
Wclosed = feedback(Wopen, 1);

popen = pole(Wopen)
pclosed = pole(Wclosed)
Dopen = [0.1 0.1 1 1];
Dclosed = [0.1 0.1 1 2];
fprintf('Open: a2*a1=%.2f, a3*a0=%.2f\n', Dopen(2)*Dopen(3), Dopen(1)*Dopen(4));
fprintf('Closed: a2*a1=%.2f, a3*a0=%.2f\n', Dclosed(2)*Dclosed(3), Dclosed(1)*Dclosed(4));

figure;
nyquist(Wopen);
grid on;
hold on;
plot(-1, 0, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
xlim([-5 5]);
ylim([-5 5]);
title('Годограф Найквиста разомкнутой системы');

D = [0.1 0.1 1 2];
w = 0:0.002:20;
Djw = polyval(D, 1i*w);
figure;
plot(real(Djw), imag(Djw), 'LineWidth', 1.8);
grid on;
title('Кривая Михайлова замкнутой системы');
xlabel('Re D(j\omega)');
ylabel('Im D(j\omega)');

figure;
step(Wopen, Wclosed);
grid on;
title('Переходные процессы');
legend('Разомкнутая система', 'Замкнутая система');
