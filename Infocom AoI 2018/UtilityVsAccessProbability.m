% Figure 3
clc;
clear all;
close all;

beta = 0.1;               
l_idle = beta;          %length of an idle slot
l_col = 1+beta;         %length of collision slot
M = [1 2 5];                % number of DSRC nodes
N = [1 2 5];              % number of WiFi nodes
t_w = 0.2;
t_d = 1e-2:1e-2:0.99;

% Generating age plots
for m = 1:numel(M)
    age = [];
    for n = 1:numel(N)
        age = [age ; ((((1-(((1-t_d).^M(m)).*((1-t_w).^N(n)))+beta)./(t_d.*((1-t_d).^(M(m)-1)).*((1-t_w).^N(n))))+(beta/2)+...
                (((1+beta)*(1-(((1-t_d).^M(m)).*((1-t_w).^N(n)))))./(2*(1-(((1-t_d).^M(m)).*((1-t_w).^N(n)))+beta)))))];
    end
    AGE{m} = age;
end

for m = 1:numel(M)
    figure;
    plot(t_d,AGE{m}(1,:),'-*r','LineWidth',2,'MarkerSize',22);
    hold on;
    plot(t_d,AGE{m}(2,:),'-+k','LineWidth',2,'MarkerSize',22);
    plot(t_d,AGE{m}(3,:),'-ob','LineWidth',2,'MarkerSize',22);
    hold off;
    set(gca,'FontSize',42,'FontWeight','bold');
    fig_name = sprintf('Fix_N_D_%d_vary_N_W.pdf',M(m));
    xlabel('\tau_D','fontweight','bold','fontsize',42);
    ylabel('AoI (\Delta)','fontweight','bold','fontsize',42);
    addpath 'legendflex'
    addpath 'setgetpos_V1.2'
    legendflex(gca,{'N_W = 1', 'N_W = 2','N_W = 5'},'ncol',3,'fontsize',40,'FontWeight','bold','anchor',{'n','n'});
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', get(0,'Screensize'));
    if m==1
        ylim([0 40]);
    elseif m==2
        ylim([0 70]);
    elseif m==3
        ylim([0 150]);
    end
    cd 'Figures'
    addpath '..\Export_fig'
    addpath '..\setgetpos_V1.2'
    export_fig(fig_name);
    cd '..\'
    close all;
end

M = [1 2 5];
N = [1 2 5];
t_w = 0.01:0.01:0.99;
t_d = 0.2;
% Generating throughput plots
for n = 1:numel(N)
    thr = [];
    for m = 1:numel(M)
        thr = [thr ; (N(n).*t_w.*((1-t_w).^(N(n)-1)).*((1-t_d)^M(m))*(1+beta))./(1-(((1-t_d).^M(m))*((1-t_w).^N(n)))+beta)];
    end
    THR{n} = thr;
end

for n = 1:numel(N)
    figure;
    plot(t_w,THR{n}(1,:),'-*r','LineWidth',2,'MarkerSize',22);
    hold on;
    plot(t_w,THR{n}(2,:),'-+k','LineWidth',2,'MarkerSize',22);
    plot(t_w,THR{n}(3,:),'-ob','LineWidth',2,'MarkerSize',22);
    hold off;
    set(gca,'FontSize',42,'FontWeight','bold');
    fig_name = sprintf('Fix_N_W_%d_vary_N_D.pdf',N(n));
    xlabel('\tau_W','fontweight','bold','fontsize',42);
    ylabel('Throughput (T)','fontweight','bold','fontsize',42);
    addpath 'legendflex'
    addpath 'setgetpos_V1.2'
    legendflex(gca,{'N_D = 1', 'N_D = 2','N_D = 5'},'ncol',3,'fontsize',40,'FontWeight','bold','anchor',{'n','n'});
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', get(0,'Screensize'));
    if n == 1
        ylim([0 1]);
    elseif n == 2
        ylim([0 1]);
    elseif n == 3
        ylim([0 1]);
    end
    cd 'Figures'
    addpath '..\Export_fig'
    addpath '..\setgetpos_V1.2'
    export_fig(fig_name);
    cd '..\'
    close all;
end