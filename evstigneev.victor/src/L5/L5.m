clc; clear; close all;

K0 = 200;
T0 = 0.2;
Mmax = 1.5;
band = 0.05;

Phi0 = tf(K0, [T0 1 K0]);  
wc0 = sqrt((-1 + sqrt(1 + 4*T0^2*K0^2))/(2*T0^2));
gamma0 = 90 - atand(T0*wc0);
fprintf('До коррекции: wc = %.4f рад/с, запас по фазе = %.2f град.\n', wc0, gamma0);
fprintf('Миноры Гурвица: %.0f и %.0f\n', 1, K0);

gamma = asind(1/Mmax);
w = 0.001:0.001:400;  
jw = 1i*w;
W0 = K0./(jw.*(1 + T0*jw));

left = 100;
right = 150;
searchIterations = 25;
for i = 1:searchIterations
    wc = (left + right)/2;
    [K,T1,T2] = corrector(wc,K0,T0,gamma);
    Wj = W0.*K.*(1 + T1*jw)./(1 + T2*jw);
    M = max(abs(Wj./(1 + Wj)));
    
    if M <= Mmax
        left = wc;
    else
        right = wc;
    end
end

wc = left;
[K,T1,T2,delta] = corrector(wc,K0,T0,gamma);
wz = 1/T1;
wp = 1/T2;
fprintf('\nВыбрано wc = %.4f рад/с, delta = %.2f град.\n', wc, delta);
fprintf('K = %.6f, T1 = %.8f с, T2 = %.8f с\n', K,T1,T2);
fprintf('wz = %.2f рад/с, wp = %.2f рад/с\n', wz,wp);
fprintf('Срез по асимптотам: %.4f и %.4f рад/с\n', sqrt(K0/T0),K0*K*T1/T0);
disp('Корректор Wку');
disp(tf(K*[T1 1], [T2 1]));

t0 = 0:0.00004:2;
t1 = 0:0.000003:0.15;
[y0,sigma0,Tpp0] = stepValues(Phi0,t0,band);
M0 = max(abs(W0./(1 + W0)));
fprintf('\n   wc, рад/с         M    sigma, %%     Tпп, с\n');
for wcTest = [75 100 wc 110]
    [Kt,T1t,T2t] = corrector(wcTest,K0,T0,gamma);
    Wtest = W0.*Kt.*(1 + T1t*jw)./(1 + T2t*jw);
    Mtest = max(abs(Wtest./(1 + Wtest)));
    num = K0*Kt*[T1t 1];
    den = [T0*T2t, T0+T2t, 1+K0*Kt*T1t, K0*Kt];
    [~,sigma,Tpp] = stepValues(tf(num,den),t1,band);
    fprintf('%12.4f %9.5f %11.4f %10.5f\n',wcTest,Mtest,sigma,Tpp);
end

num = K0*K*[T1 1];
den = [T0*T2, T0+T2, 1+K0*K*T1, K0*K];
Phi1 = tf(num,den);
disp('Замкнутая система после коррекции');
disp(Phi1);
D2 = den(2)*den(3) - den(1)*den(4);
fprintf('Миноры Гурвица: %.5f, %.5f, %.2f\n',den(2),D2,den(4)*D2);

[y1,sigma1,Tpp1] = stepValues(Phi1,t1,band);
Wj = W0.*K.*(1 + T1*jw)./(1 + T2*jw);
M1 = max(abs(Wj./(1 + Wj)));
fprintf('До:    M = %.5f, sigma = %.4f %%, Tпп = %.5f с\n',M0,sigma0,Tpp0);
fprintf('После: M = %.5f, sigma = %.4f %%, Tпп = %.5f с\n',M1,sigma1,Tpp1);

w = logspace(-0.5,3.5,2000);
jw = 1i*w;
W0 = K0./(jw.*(1 + T0*jw));
C = K*(1 + T1*jw)./(1 + T2*jw);
Wj = W0.*C;
L0as = 20*log10(K0./w) - 20*log10(max(1,T0*w));
Ljas = L0as + 20*log10(K) + 20*log10(max(1,T1*w)) - 20*log10(max(1,T2*w));

figure('Name','Разомкнутые контуры');
subplot(2,1,1);
semilogx(w,20*log10(abs(W0)),w,20*log10(abs(Wj)),w,L0as,'--',w,Ljas,'--');
grid on; ylabel('L, дБ'); xlabel('Частота, рад/с');
legend('Исходная','Желаемая','Асимптота исходной','Асимптота желаемой');
subplot(2,1,2);
semilogx(w,angle(W0)*180/pi,w,angle(Wj)*180/pi);
grid on; ylabel('Фаза, град.'); xlabel('Частота, рад/с');
legend('Исходная','Желаемая');

figure('Name','ЛАЧХ корректора');
semilogx(w,Ljas-L0as,'--',w,20*log10(abs(C)));
grid on; ylabel('L, дБ'); xlabel('Частота, рад/с');
legend('Разность асимптот','Точная ЛАЧХ');

F0 = W0./(1 + W0);
F1 = Wj./(1 + Wj);
figure('Name','Частотные характеристики замкнутых систем');
subplot(2,1,1);
semilogx(w,20*log10(abs(F0)),w,20*log10(abs(F1)));
grid on; ylabel('L, дБ'); xlabel('Частота, рад/с');
legend('До коррекции','После коррекции');
subplot(2,1,2);
semilogx(w,angle(F0)*180/pi,w,angle(F1)*180/pi);
grid on; ylabel('Фаза, град.'); xlabel('Частота, рад/с');

figure('Name','Переходные характеристики');
subplot(1,2,1);
plot(t0,y0); grid on; title('До коррекции');
yline(1-band,':'); yline(1+band,':'); xline(Tpp0,':');
xlabel('Время, с'); ylabel('h(t)');
subplot(1,2,2);
plot(t1,y1); grid on; title('После коррекции');
yline(1-band,':'); yline(1+band,':'); xline(Tpp1,':');
xlabel('Время, с'); ylabel('h(t)');

function [K,T1,T2,delta] = corrector(wc,K0,T0,gamma)
    delta = gamma - (90 - atand(T0*wc));
    alpha = (1+sind(delta))/(1-sind(delta));
    T2 = 1/(wc*sqrt(alpha));
    T1 = alpha*T2;
    K = wc*sqrt(1+(T0*wc)^2)/(K0*sqrt(alpha));
end

function [y,sigma,Tpp] = stepValues(Phi,t,band)
    y = step(Phi,t);
    sigma = 100*(max(y) - 1);
    e = abs(y - 1) - band;
    i = find(e > 0,1,'last');
    Tpp = t(i) + (t(i+1)-t(i))*e(i)/(e(i)-e(i+1));
end
