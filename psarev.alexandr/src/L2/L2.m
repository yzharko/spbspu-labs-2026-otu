clc;
clear;
close all;

syms s

W1 = 1 / (s + 1);
W2 = 2;
W3 = 5 / (s^2 + 3*s + 2);
W4 = 4;

W_final = (W1*(1 + W3) - W3*W4) / (1 + W3 + W3*W4*(1 + W2 + W1*W2));

W_simpl = simplify(W_final)

[num, den] = numden(W_simpl);
num_coeffs = sym2poly(num)
den_coeffs = sym2poly(den)

names = {'скачок', 'импульс', 'синусоида'};
if exist('scheme_L2.slx', 'file')
    load_system('scheme_L2');
    for k = 1:3
        set_param('scheme_L2/u_type', 'Value', num2str(k));
        out = sim('scheme_L2');
        y1 = out.logsout.get('y').Values;
        y2 = out.logsout.get('y_eq').Values;
        figure;
        plot(y1.Time, y1.Data, 'b', y2.Time, y2.Data, 'r--', 'LineWidth', 1.5);
        grid on; legend('Исходная схема', 'W\_check');
        xlabel('t, c'); ylabel('y'); title([names{k}]);
        fprintf('%s: max |y - y_eq| = %g\n', names{k}, ...
            max(abs(y1.Data - interp1(y2.Time, y2.Data, y1.Time))));
    end
    set_param('scheme_L2/u_type', 'Value', '1');
end
