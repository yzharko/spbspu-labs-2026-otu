clear; clc; close all;

s = tf('s');
K = 200;
Tp = [0.1 0.02 0.01];
M_req = 1.35;
wc = 12;

W = K/(s*(0.1*s+1)*(0.02*s+1)*(0.01*s+1));
wp = 1./Tp;

q = M_req/(M_req-1);
w2 = wc/q;
w1 = wc*w2/K;

wp_cancel = wp(wp <= wc);
wp_keep = wp(wp > wc);
wp_extra = repmat(wp(end), 1, numel(wp_cancel));

Wku = (s/w2+1)/(s/w1+1);
for p = wp_cancel, Wku = Wku*(s/p+1); end
for p = wp_extra,  Wku = Wku/(s/p+1); end
Wku = minreal(Wku);

Wsk = minreal(W*Wku);
T0 = feedback(W,1);
T = feedback(Wsk,1);

w = logspace(-1,3,3000);
Lc = @(c) sum(20*log10(max(w(:)./c(:).',1)),2).';
Lsrc = 20*log10(K./w) - Lc(wp);
Lzh = 20*log10(K./w) + Lc(w2) - Lc([w1 wp_keep wp_extra]);
Lku = Lzh - Lsrc;

f1 = figure;
semilogx(w,Lsrc,'b',w,Lzh,'r',w,Lku,'g--','LineWidth',1.5);
grid on; yline(0,'k'); xline(wc,':k','\omega_c');
legend('L_{исх}','L_{ж}','L_{ку}','Location','southwest');
xlabel('\omega, рад/с'); ylabel('L, дБ'); title('Асимптотические ЛАЧХ');

f2 = figure;
step(T0,0:0.001:0.5); grid on;
title('Переходная характеристика исходной САУ');

f3 = figure;
step(T,0:0.001:2); grid on;
title('Переходная характеристика после коррекции');

f4 = figure;
margin(Wsk); grid on;

f5 = figure;
bode(W,'b',Wsk,'r',Wku,'g--'); grid on;
legend('W','W\cdot W_{ку}','W_{ку}');

saveas(f1,'fig_lachx.png');
saveas(f2,'fig_step_before.png');
saveas(f3,'fig_step_after.png');
saveas(f4,'fig_margin.png');
saveas(f5,'fig_bode.png');
