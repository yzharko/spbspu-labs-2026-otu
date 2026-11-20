clc;
clear;
close all;

s = tf('s');

D_open = [2 3 3 0.5 0.5 1];
D_closed = [2 3 3 0.5 0.5 11];

W_open = tf(10, D_open);
W_closed = feedback(W_open, 1);

polynomials = {D_open, D_closed};
PZ = zeros(1, 2);

for m = 1:2
    a = polynomials{m};
    n = length(a) - 1;
    H = zeros(n);

    for row = 1:n
        for col = 1:n
            idx = 2 * col - row + 1;

            if idx >= 1 && idx <= length(a)
                H(row, col) = a(idx);
            end
        end
    end

    minors = zeros(1, n);

    for k = 1:n
        minors(k) = det(H(1:k, 1:k));
    end

    poles = roots(a);
    PZ(m) = sum(real(poles) > 1e-8);

    fprintf('\nSystem %d\n', m);
    disp(H);
    disp(minors);
    disp(poles);

    fprintf('Hurwitz stable: %d\n', all(minors > 0));
end

w = [0, logspace(-5, 4, 6000)];

[re, im] = nyquist(W_open, w);
re = squeeze(re);
im = squeeze(im);

figure('Color', 'w');
plot(re, im, 'b', re, -im, 'b');
hold on;
plot(-1, 0, 'r+');
grid on;

xlabel('Re W(j\omega)');
ylabel('Im W(j\omega)');
title('Nyquist');

figure('Color', 'w');
plot(re, im, 'b', re, -im, 'b');
hold on;
plot(-1, 0, 'r+');
grid on;

xlim([-2 2]);
ylim([-2 2]);
xlabel('Re W(j\omega)');
ylabel('Im W(j\omega)');
title('Nyquist near (-1, 0)');

wp = logspace(-6, 5, 20000);
wf = [-fliplr(wp), 0, wp];

phase = unwrap(angle(1 + 10 ./ polyval(D_open, 1i * wf)));
N = round((phase(end) - phase(1)) / (2 * pi));

fprintf('Clockwise N=%d, P=%d, Z=%d\n', N, PZ(1), PZ(1) + N);

for maxw = [3, 1.3]
    w = linspace(0, maxw, 5000);
    D = polyval(D_closed, 1i * w);

    figure('Color', 'w');
    plot(real(D), imag(D), 'LineWidth', 1.5);
    grid on;

    xlabel('Re D(j\omega)');
    ylabel('Im D(j\omega)');
    title('Mikhailov');
end

w = [0, logspace(-6, 5, 20000)];
ph = unwrap(angle(polyval(D_closed, 1i * w)));

fprintf('Mikhailov phase change / pi = %.6f\n', ...
    (ph(end) - ph(1)) / pi);

figure('Color', 'w');

subplot(2, 1, 1);
step(W_open, 0:0.01:24);
grid on;
title('Open loop');

subplot(2, 1, 2);
step(W_closed, 0:0.01:12);
grid on;
title('Closed loop');
