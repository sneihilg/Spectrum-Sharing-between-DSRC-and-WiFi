% plotting best response curve for 2 DSRC networks without cost 
clc;
clear all;
close all;

beta = 0.001;     
l_idle = beta;          %length of an idle slot
l_col = 1+beta;         %length of collision slot
M = [5 5 5];                  % number of WiFi nodes in network 1
N = [1 2 5];                  % number of WiFi nodes in network 2
tau_w1 = 1e-2:1e-2:0.99;
tau_w2 = 1e-2:1e-2:0.99;

options = optimset('Algorithm','interior-point','TolX',1e-10,...
     'TolFun',1e-8,'TolCon',1e-10,'MaxFunEval',1e6,'MaxIter',1e6,'Display','Iter','InitBarrierParam',1e-2);

optimal_tau_w1 = [];
optimal_tau_w2 = [];
NE = [];
for m = 1:numel(M)
    optim_tau_w1 = [];
    optim_thr1 = [];
    for i = 1:numel(tau_w2)
        t_w_o = 1e-2;
        lb = 1e-2;
        ub = 0.99;
        [t_w1,fval,exitflag] = fmincon(@(t_w1)optimizing_net_thr(beta,t_w1,tau_w2(i),N(m),M(m)),t_w_o,[],[],[],[],lb,ub,[],options);
        optim_tau_w1 = [optim_tau_w1 t_w1];
        optim_thr1 = [optim_thr1 -fval];
    end
    optimal_tau_w1 = [optimal_tau_w1;optim_tau_w1];
    
    optim_tau_w2 = [];
    optim_thr2 = [];
    for i = 1:numel(tau_w1)
        t_w_o = 1e-2;
        lb = 1e-2;
        ub = 0.99;
        [t_w2,fval,exitflag] = fmincon(@(t_w2)optimizing_net_thr(beta,t_w2,tau_w1(i),M(m),N(m)),t_w_o,[],[],[],[],lb,ub,[],options);
        optim_tau_w2 = [optim_tau_w2 t_w2];
        optim_thr2 = [optim_thr2 -fval];
    end
    optimal_tau_w2 = [optimal_tau_w2;optim_tau_w2];
    % find point of intersection of plots
    [xout,yout] = intersections(optimal_tau_w1(m,:),tau_w2,tau_w1,optimal_tau_w2(m,:),1);
    NE{m} = [xout yout];
    
    figure;
    plot(optimal_tau_w1(m,:),tau_w2,'-r','LineWidth',8,'MarkerSize',18);
    hold on;
    plot(tau_w1,optimal_tau_w2(m,:),':k','LineWidth',8,'MarkerSize',18);
    plot(xout,yout,'ob','LineWidth',8,'MarkerSize',18);
    for i = 1:numel(xout)
        if M(m)==1 || N(m)==1
            text(xout(i),yout(i),['(' num2str(xout(i),'%.2f') ',' num2str(yout(i),'%.2f') ')'],'FontSize', 42,...
                'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top');
        else
            text(xout(i),yout(i),['(' num2str(xout(i),'%.2f') ',' num2str(yout(i),'%.2f') ')'],'FontSize', 42,...
            'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','bottom');
        end
    end
    hold off;
    set(gca,'FontSize',42,'FontWeight','bold');
    fig_name = sprintf('M_%d_N_%d_WiFi.pdf',M(m),N(m));
%     title(sprintf('N_D = %d, N_W = %d',M(m),N(m)),'fontsize',36);
    xlabel('\tau_W_1','fontweight','bold','fontsize',42);
    ylabel('\tau_W_2','fontweight','bold','fontsize',42);
    if M(m)<5
        legend('WiFi Network 1','WiFi Network 2','Location','northwest');
    else
        legend('WiFi Network 1','WiFi Network 2','Location','northeast');
    end
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', get(0,'Screensize'));
    cd 'Figures'
    addpath '..\Export_fig'
    export_fig(fig_name);
    cd '..\'
    close all;
end

for m=1:numel(M)
    t_w1 = NE{m}(:,1);
    t_w2 = NE{m}(:,2);
    thr_1 = (t_w1.*((1-t_w1).^(M(m)-1)).*((1-t_w2).^N(m)).*(1+beta))./(1-(((1-t_w1).^M(m)).*((1-t_w2).^N(m)))+beta);
    thr_2 = (t_w2.*((1-t_w2).^(N(m)-1)).*((1-t_w1).^M(m)).*(1+beta))./(1-(((1-t_w2).^N(m)).*((1-t_w1).^M(m)))+beta);
    disp('Equilibrium points');
    disp([t_w1 t_w2]);
    disp('Corresponding throughputs');
    disp([thr_1 thr_2]);
end
