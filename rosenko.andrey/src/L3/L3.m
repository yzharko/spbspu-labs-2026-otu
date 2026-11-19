clc;
clear;
close all;
set(0,'DefaultFigureVisible','off');

K0 = 0.8;
T0 = 0.5;
T1 = 0.7;
T2 = 0.5;
T1a = 1.25;
K_values = [0.4 0.8 1.6];
T_values = [0.1 0.25 0.5 1.0 2.5];
w = logspace(-2,3,900);
s = tf('s');

names = {'Integrator','Ideal differentiator','Gain','First-order lag', ...
    'Second-order aperiodic','Oscillatory','Conservative', ...
    'Lag-integrator','Real differentiator','First-order lead','Isodromic'};

outDir = 'L3_results';
if ~exist(outDir,'dir')
    mkdir(outDir);
end

% Base frequency responses for all 11 links.
for i = 1:11
    W = localW(i,K0,T0,T1,T2,T1a,s);
    [mag,phase] = bode(W,w);
    mag = squeeze(mag);
    phase = squeeze(phase);
    [re,im] = nyquist(W,w);
    re = squeeze(re); im = squeeze(im);

    f = figure('Color','w','Position',[100 100 1100 420]);
    subplot(1,3,1);
    semilogx(w,20*log10(mag),'LineWidth',1.3); grid on;
    xlabel('w, rad/s'); ylabel('L(w), dB'); title([names{i} ' magnitude']);
    subplot(1,3,2);
    semilogx(w,phase,'LineWidth',1.3); grid on;
    xlabel('w, rad/s'); ylabel('Phase, deg'); title([names{i} ' phase']);
    subplot(1,3,3);
    plot(re,im,'LineWidth',1.3); hold on; plot(re,-im,'--','LineWidth',1.0);
    plot(0,0,'k.'); grid on; axis tight;
    xlabel('Real'); ylabel('Imaginary'); title([names{i} ' Nyquist']);
    print(f,fullfile(outDir,sprintf('base_%02d.png',i)),'-dpng','-r150');
    close(f);
end

% Gain variation for each link.
for i = 1:11
    f = figure('Color','w','Position',[100 100 1100 700]);
    ax1 = subplot(3,1,1); hold(ax1,'on'); grid(ax1,'on');
    ax2 = subplot(3,1,2); hold(ax2,'on'); grid(ax2,'on');
    ax3 = subplot(3,1,3); hold(ax3,'on'); grid(ax3,'on');
    for k = K_values
        W = localW(i,k,T0,T1,T2,T1a,s);
        [mag,phase] = bode(W,w); mag=squeeze(mag); phase=squeeze(phase);
        [re,im] = nyquist(W,w); re=squeeze(re); im=squeeze(im);
        semilogx(ax1,w,20*log10(mag),'LineWidth',1.2);
        semilogx(ax2,w,phase,'LineWidth',1.2);
        plot(ax3,re,im,'LineWidth',1.2); plot(ax3,re,-im,'--','LineWidth',0.8);
    end
    xlabel(ax1,'w, rad/s'); ylabel(ax1,'L(w), dB'); title(ax1,[names{i} ' - gain variation']);
    xlabel(ax2,'w, rad/s'); ylabel(ax2,'Phase, deg');
    xlabel(ax3,'Real'); ylabel(ax3,'Imaginary'); title(ax3,'Nyquist');
    legend(ax1,{'K=0.4','K=0.8','K=1.6'},'Location','best');
    print(f,fullfile(outDir,sprintf('gain_%02d.png',i)),'-dpng','-r150');
    close(f);
end

% Common time-constant variation for links 4, 7, 8, 9, 10 and 11.
indices = [4 7 8 9 10 11];
for q = 1:length(indices)
    i = indices(q);
    f = figure('Color','w','Position',[100 100 1100 700]);
    ax1 = subplot(3,1,1); hold(ax1,'on'); grid(ax1,'on');
    ax2 = subplot(3,1,2); hold(ax2,'on'); grid(ax2,'on');
    ax3 = subplot(3,1,3); hold(ax3,'on'); grid(ax3,'on');
    for T = T_values
        W = localW(i,K0,T,T1,T2,T1a,s);
        [mag,phase] = bode(W,w); mag=squeeze(mag); phase=squeeze(phase);
        [re,im] = nyquist(W,w); re=squeeze(re); im=squeeze(im);
        semilogx(ax1,w,20*log10(mag),'LineWidth',1.1);
        semilogx(ax2,w,phase,'LineWidth',1.1);
        plot(ax3,re,im,'LineWidth',1.1); plot(ax3,re,-im,'--','LineWidth',0.8);
    end
    xlabel(ax1,'w, rad/s'); ylabel(ax1,'L(w), dB'); title(ax1,[names{i} ' - T variation']);
    xlabel(ax2,'w, rad/s'); ylabel(ax2,'Phase, deg');
    xlabel(ax3,'Real'); ylabel(ax3,'Imaginary'); title(ax3,'Nyquist');
    legend(ax1,{'T=0.1','T=0.25','T=0.5','T=1','T=2.5'},'Location','best');
    print(f,fullfile(outDir,sprintf('time_%02d.png',i)),'-dpng','-r150');
    close(f);
end

% Second-order elements: vary T1 (or T1a) and T2 separately.
% The selected values preserve the required aperiodic/oscillatory conditions.
second = [5 6];
T1a_values = [1.0 1.25 1.5 2.0 2.5];
T2_aper_values = [0.25 0.4 0.5 0.6];
T1_osc_values = [0.1 0.3 0.5 0.7];
T2_osc_values = [0.4 0.5 0.7 1.0];

for q = 1:length(second)
    i = second(q);
    f = figure('Color','w','Position',[100 100 1100 700]);
    ax1 = subplot(2,3,1); hold(ax1,'on'); grid(ax1,'on');
    ax2 = subplot(2,3,2); hold(ax2,'on'); grid(ax2,'on');
    ax3 = subplot(2,3,3); hold(ax3,'on'); grid(ax3,'on');
    ax4 = subplot(2,3,4); hold(ax4,'on'); grid(ax4,'on');
    ax5 = subplot(2,3,5); hold(ax5,'on'); grid(ax5,'on');
    ax6 = subplot(2,3,6); hold(ax6,'on'); grid(ax6,'on');
    if i == 5
        values = T1a_values;
        for v = values
            W = localW(5,K0,T0,T1,T2,v,s);
            [mag,phase] = bode(W,w); mag=squeeze(mag); phase=squeeze(phase);
            [re,im] = nyquist(W,w); re=squeeze(re); im=squeeze(im);
            semilogx(ax1,w,20*log10(mag),'LineWidth',1.1);
            semilogx(ax2,w,phase,'LineWidth',1.1);
            plot(ax3,re,im,'LineWidth',1.1); plot(ax3,re,-im,'--','LineWidth',0.8);
        end
        values = T2_aper_values;
        for v = values
            W = localW(5,K0,T0,T1,v,T1a,s);
            [mag,phase] = bode(W,w); mag=squeeze(mag); phase=squeeze(phase);
            [re,im] = nyquist(W,w); re=squeeze(re); im=squeeze(im);
            semilogx(ax4,w,20*log10(mag),'LineWidth',1.1);
            semilogx(ax5,w,phase,'LineWidth',1.1);
            plot(ax6,re,im,'LineWidth',1.1); plot(ax6,re,-im,'--','LineWidth',0.8);
        end
        legend(ax1,{'T1a=1.0','T1a=1.25','T1a=1.5','T1a=2.0','T1a=2.5'},'Location','best');
        legend(ax4,{'T2=0.25','T2=0.4','T2=0.5','T2=0.6'},'Location','best');
    else
        values = T1_osc_values;
        for v = values
            W = localW(6,K0,T0,v,T2,T1a,s);
            [mag,phase] = bode(W,w); mag=squeeze(mag); phase=squeeze(phase);
            [re,im] = nyquist(W,w); re=squeeze(re); im=squeeze(im);
            semilogx(ax1,w,20*log10(mag),'LineWidth',1.1);
            semilogx(ax2,w,phase,'LineWidth',1.1);
            plot(ax3,re,im,'LineWidth',1.1); plot(ax3,re,-im,'--','LineWidth',0.8);
        end
        values = T2_osc_values;
        for v = values
            W = localW(6,K0,T0,T1,v,T1a,s);
            [mag,phase] = bode(W,w); mag=squeeze(mag); phase=squeeze(phase);
            [re,im] = nyquist(W,w); re=squeeze(re); im=squeeze(im);
            semilogx(ax4,w,20*log10(mag),'LineWidth',1.1);
            semilogx(ax5,w,phase,'LineWidth',1.1);
            plot(ax6,re,im,'LineWidth',1.1); plot(ax6,re,-im,'--','LineWidth',0.8);
        end
        legend(ax1,{'T1=0.1','T1=0.3','T1=0.5','T1=0.7'},'Location','best');
        legend(ax4,{'T2=0.4','T2=0.5','T2=0.7','T2=1.0'},'Location','best');
    end
    xlabel(ax1,'w, rad/s'); ylabel(ax1,'L(w), dB'); title(ax1,'Coefficient 1: magnitude');
    xlabel(ax2,'w, rad/s'); ylabel(ax2,'Phase, deg'); title(ax2,'Coefficient 1: phase');
    xlabel(ax3,'Real'); ylabel(ax3,'Imaginary'); title(ax3,'Coefficient 1: Nyquist');
    xlabel(ax4,'w, rad/s'); ylabel(ax4,'L(w), dB'); title(ax4,'Coefficient 2: magnitude');
    xlabel(ax5,'w, rad/s'); ylabel(ax5,'Phase, deg'); title(ax5,'Coefficient 2: phase');
    xlabel(ax6,'Real'); ylabel(ax6,'Imaginary'); title(ax6,'Coefficient 2: Nyquist');
    print(f,fullfile(outDir,sprintf('second_order_%02d.png',i)),'-dpng','-r150');
    close(f);
end

fprintf('Base parameters: K=%.2f, T=%.2f, T1=%.2f, T2=%.2f, T1a=%.2f\n',K0,T0,T1,T2,T1a);
fprintf('K values: 0.4, 0.8, 1.6\n');
fprintf('Frequency plots and parameter studies were saved to %s\n',outDir);
pngFiles = dir(fullfile(outDir,'*.png'));
if ~isempty(pngFiles)
    fileList = cell(1,numel(pngFiles));
    for q = 1:numel(pngFiles)
        fileList{q} = fullfile(outDir,pngFiles(q).name);
    end
    zipFile = 'L3_graphs.zip';
    if exist(zipFile,'file'), delete(zipFile); end
    zip(zipFile,fileList);
end
set(0,'DefaultFigureVisible','on');

function W = localW(i,Kx,Tx,T1x,T2x,T1ax,s)
    if i == 1
        W = Kx/s;
    elseif i == 2
        W = Kx*s;
    elseif i == 3
        W = tf(Kx);
    elseif i == 4
        W = Kx/(Tx*s+1);
    elseif i == 5
        W = Kx/(T2x^2*s^2+T1ax*s+1);
    elseif i == 6
        W = Kx/(T2x^2*s^2+T1x*s+1);
    elseif i == 7
        W = Kx/(Tx*s^2+1);
    elseif i == 8
        W = Kx/(s*(Tx*s+1));
    elseif i == 9
        W = Kx*s/(Tx*s+1);
    elseif i == 10
        W = Kx*(Tx*s+1);
    elseif i == 11
        W = Kx*(Tx*s+1)/s;
    else
        error('Unknown link index');
    end
end
