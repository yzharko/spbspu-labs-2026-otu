clear; clc; close all;

syms s t

K = sym("1.7");
T = sym("0.3");
T1 = sym("0.4");
T2 = sym("0.6");

Ws = {
    K / s
    K * s
    K
    K / (T * s + 1)
    K / (T2^2 * s^2 + T1 * s + 1)
    K / (T * s^2 + 1)
    K / (s * (T * s + 1))
    K * s / (T * s + 1)
    K * (T * s + 1)
    K * (T * s + 1) / s
};

U_step = laplace(heaviside(t), t, s);
U_dirac = laplace(dirac(t), t, s);

t_range = 0:0.01:10;

plot_all("Ступенька", Ws, U_step, s, t, t_range, 5, 2);
plot_all("Дирак", Ws, U_dirac, s, t, t_range, 5, 2);

% Теперь задание

% 1 Варьировать буду осцилятор.
coeffs = {
    K * 2 / (T2^2 * s^2 + T1 * s + 1)
    K / (T2^2 * s^2 + T1 * s + 1)
    K / 2 / (T2^2 * s^2 + T1 * s + 1)
};
plot_all("Коэффициенты", coeffs, U_dirac, s, t, t_range, 3, 1);

% 2
aperiodic = {
    K / (T * 5 * s + 1)
    K / (T * 2 * s + 1)
    K / (T * s + 1)
    K / (T / 2 * s + 1)
    K / (T / 5 * s + 1)
};
plot_all("Апериодическое", aperiodic, U_step, s, t, t_range, 5, 1);

% 3
oscilator = {
    K / (T2^2 * s^2 + T1 * s + 1)
    K / (T2^2 * s^2 + T1 * 0.5 * s + 1)
    K / (T2^2 * s^2 + T1 * 0.2 * s + 1)
    K / (T2^2 * s^2 + T1 * 0.1 * s + 1)
    K / (T2^2 * s^2 + T1 * 0.0 * s + 1)
};
plot_all("Осцилятор", oscilator, U_step, s, t, t_range, 5, 1);

function plot_all(name, Ws, U, s, t, t_range, r, c)
    figure("Name", name);
    for i = 1:numel(Ws)
        F = Ws{i} * U;
        y = matlabFunction(ilaplace(F, s, t), "Vars", t);
        y_range = y(t_range);
        subplot(r, c, i);
        plot(t_range, y_range);
        grid on;
        title(string(F));
    end
end
