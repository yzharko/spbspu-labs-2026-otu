clear; clc; close all;

syms u y e1 e2 e3 y1 y2 y4 W1 W2 W3 W4 s

%{
e1 = u - y4
y1 = W1 * e1
e2 = y1 - y
y2 = W2 * e2
e3 = y2 + e1
y = W3 * e3
y4 = W4 * y

W1 = 1 / (s + 1)
W2 = 2
W3 = 5 / (s^2 + 3 * s + 2)
W4 = 4
%}

steps = {};
steps{1} = u - y4 == e1;
steps{2} = subs(steps{1}, e1, y1 / W1);
steps{3} = subs(steps{2}, y1, e2 + y);
steps{4} = subs(steps{3}, e2, y2 / W2);
steps{5} = subs(steps{4}, y2, e3 - e1);
steps{6} = subs(steps{5}, e3, y / W3);
steps{7} = subs(steps{6}, e1, u - y4);
steps{8} = subs(steps{7}, y4, W4 * y);
steps{9} = solve(steps{8}, y) / u;
steps{10} = simplifyFraction(steps{9});
steps{11} = subs(steps{10}, W1, 1 / (s + 1));
steps{12} = subs(steps{11}, W2, 2);
steps{13} = subs(steps{12}, W3, 5 / (s^2 + 3 * s + 2));
steps{14} = subs(steps{13}, W4, 4);
steps{15} = simplifyFraction(steps{14}, "Expand", true);

for i = 1:numel(steps)
    disp(steps{i})
end

% num: 5, 15
% den: 1, 4, 35, 72