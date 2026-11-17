clear;
clc;
close all;

%% Передаточные функции звеньев

W1 = tf(1, [1 1]);
W2 = tf(2);
W3 = tf(5, [1 3 2]);
W4 = tf(4);

%% Соединения исходной структурной схемы

W1.InputName = 'x';
W1.OutputName = 'afterW1';

W2.InputName = 'y';
W2.OutputName = 'afterW2';

W3.InputName = 'beforeW3';
W3.OutputName = 'afterW3';

W4.InputName = 'afterW3';
W4.OutputName = 'feedbackSignal';

firstSum = sumblk('x = u - feedbackSignal');
secondSum = sumblk('y = afterW1 - feedbackSignal');
thirdSum = sumblk('beforeW3 = x + afterW2');

originalSystem = connect(W1, W2, W3, W4, firstSum, secondSum, thirdSum, 'u', 'y');

%% Эквивалентная передаточная функция

equivalentSystem = tf([1 -17 -18], [1 4 65 102]);

disp('Исходная система:');
disp(minreal(originalSystem));

disp('Эквивалентная система:');
disp(equivalentSystem);

disp('Разность передаточных функций:');
disp(minreal(originalSystem - equivalentSystem));

%% Временной интервал

time = 0:0.01:30;

%% Сравнение переходных характеристик

originalStep = squeeze(step(originalSystem, time));
equivalentStep = squeeze(step(equivalentSystem, time));

figure('Name', 'Ступенчатое воздействие', 'NumberTitle', 'off');

plot(time, originalStep, 'LineWidth', 2);
hold on;
plot(time, equivalentStep, '--', 'LineWidth', 1.5);
grid on;

title('Сравнение переходных характеристик');
xlabel('Время t, с');
ylabel('y(t)');
legend('Исходная схема', 'Эквивалентная схема');
hold off;

%% Сравнение импульсных характеристик

originalImpulse = squeeze(impulse(originalSystem, time));
equivalentImpulse = squeeze(impulse(equivalentSystem, time));

figure('Name', 'Импульсное воздействие', 'NumberTitle', 'off');

plot(time, originalImpulse, 'LineWidth', 2);
hold on;
plot(time, equivalentImpulse, '--', 'LineWidth', 1.5);
grid on;

title('Сравнение импульсных характеристик');
xlabel('Время t, с');
ylabel('y(t)');
legend('Исходная схема', 'Эквивалентная схема');
hold off;

%% Сравнение при синусоидальном воздействии

inputSignal = sin(time);

originalSine = lsim(originalSystem, inputSignal, time);
equivalentSine = lsim(equivalentSystem, inputSignal, time);

figure('Name', 'Синусоидальное воздействие', 'NumberTitle', 'off');

plot(time, originalSine, 'LineWidth', 2);
hold on;
plot(time, equivalentSine, '--', 'LineWidth', 1.5);
grid on;

title('Сравнение реакций на синусоидальное воздействие');
xlabel('Время t, с');
ylabel('y(t)');
legend('Исходная схема', 'Эквивалентная схема');
hold off;

%% Проверка совпадения результатов

stepError = max(abs(originalStep(:) - equivalentStep(:)));
impulseError = max(abs(originalImpulse(:) - equivalentImpulse(:)));
sineError = max(abs(originalSine(:) - equivalentSine(:)));

fprintf('\nМаксимальная разность результатов:\n');
fprintf('Ступенчатое воздействие: %.3e\n', stepError);
fprintf('Импульсное воздействие: %.3e\n', impulseError);
fprintf('Синусоидальное воздействие: %.3e\n', sineError);
