% plotting best response curve for 2 DSRC networks without cost 
clc;
clear all;
close all;

beta = 0.001;     
l_idle = beta;          %length of an idle slot
l_col = 1+beta;         %length of collision slot
M = [5 5 5];                  % number of DSRC nodes in network 1
N = [1 2 5];                  % number of DSRC nodes in network 2
tau_d1 = 1e-2:1e-2:0.99;
tau_d2 = 1e-2:1e-2:0.99;

options = optimset('Algorithm','interior-point','TolX',1e-10,...
     'TolFun',1e-8,'TolCon',1e-10,'MaxFunEval',1e6,'MaxIter',1e6,'Display','Iter','InitBarrierParam',1e-12);

optimal_tau_d1 = [];
optimal_tau_d2 = [];
NE = [];
for m = 1:numel(M)
    optim_tau_d1 = [];
    optim_age_d1 = [];
    for i = 1:numel(tau_d2)
        t_d_o = 1e-2;
        lb = 1e-2;
        ub = 0.99;
        t_d1_temp = 1e-2:1e-2:0.99;
        t_d2_temp = tau_d2(i);
        age_temp = ((((1-(((1-t_d1_temp).^M(m)).*((1-t_d2_temp).^N(m)))+beta)./(t_d1_temp.*((1-t_d1_temp).^(M(m)-1)).*((1-t_d2_temp).^N(m))))+(beta/2)+...
                (((1+beta)*(1-(((1-t_d1_temp).^M(m)).*((1-t_d2_temp).^N(m)))))./(2*(1-(((1-t_d1_temp).^M(m)).*((1-t_d2_temp).^N(m)))+beta)))));
        minimum = min(age_temp);
        maximum = max(age_temp);
        [t_d1,fval,exitflag] = fmincon(@(t_d1)optimizing_net_age(beta,t_d1,tau_d2(i),minimum,maximum,M(m),N(m)),t_d_o,[],[],[],[],lb,ub,[],options);
        optim_tau_d1 = [optim_tau_d1 t_d1];
        optim_age_d1 = [optim_age_d1 fval];
    end
    optimal_tau_d1 = [optimal_tau_d1 ; optim_tau_d1];
    optim_tau_d2 = [];
    optim_age_d2 = [];
    for i = 1:numel(tau_d1)
        t_d_o = 1e-2;
        lb = 1e-2;
        ub = 0.99;
        t_d2_temp = 1e-2:1e-2:0.99;
        t_d1_temp = tau_d1(i);
        age_temp = ((((1-(((1-t_d2_temp).^N(m)).*((1-t_d1_temp).^M(m)))+beta)./(t_d2_temp.*((1-t_d2_temp).^(N(m)-1)).*((1-t_d1_temp).^M(m))))+(beta/2)+...
                (((1+beta)*(1-(((1-t_d2_temp).^N(m)).*((1-t_d1_temp).^M(m)))))./(2*(1-(((1-t_d2_temp).^N(m)).*((1-t_d1_temp).^M(m)))+beta)))));
        minimum = min(age_temp);
        maximum = max(age_temp);
        [t_d2,fval,exitflag] = fmincon(@(t_d2)optimizing_net_age(beta,t_d2,tau_d1(i),minimum,maximum,N(m),M(m)),t_d_o,[],[],[],[],lb,ub,[],options);
        optim_tau_d2 = [optim_tau_d2 t_d2];
        optim_age_d2 = [optim_age_d2 fval];
    end
    optimal_tau_d2 = [optimal_tau_d2 ; optim_tau_d2];
    
    % find point of intersection of plots
    [xout,yout] = intersections(optimal_tau_d1(m,:),tau_d2,tau_d1,optimal_tau_d2(m,:),1);
    NE{m} = [xout yout];
    figure;
    plot(optimal_tau_d1(m,:),tau_d2,'-r','LineWidth',8,'MarkerSize',18);
    hold on;
    plot(tau_d1,optimal_tau_d2(m,:),':k','LineWidth',8,'MarkerSize',18);
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
    fig_name = sprintf('M_%d_N_%d_DSRC.pdf',M(m),N(m));
    title(sprintf('N_D = %d, N_W = %d',M(m),N(m)),'fontsize',36);
    xlabel('\tau_D_1','fontweight','bold','fontsize',42);
    ylabel('\tau_D_2','fontweight','bold','fontsize',42);
    if M(m)<5
        legend('DSRC Network 1','DSRC Network 2','Location','northwest');
    else
        legend('DSRC Network 1','DSRC Network 2','Location','northeast');
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
    t_d1 = NE{m}(:,1);
    t_d2 = NE{m}(:,2);
    age_1 = ((1-(((1-t_d1).^M(m)).*((1-t_d2).^N(m)))+beta)./(t_d1.*((1-t_d1).^(M(m)-1)).*((1-t_d2).^N(m))))+(beta/2)+...
        (((1+beta).*(1-(((1-t_d1).^M(m)).*((1-t_d2).^N(m)))))./(2.*(1-(((1-t_d1).^M(m)).*((1-t_d2).^N(m)))+beta)));
    age_2 = ((1-(((1-t_d2).^N(m)).*((1-t_d1).^M(m)))+beta)./(t_d2.*((1-t_d2).^(N(m)-1)).*((1-t_d1).^M(m))))+(beta/2)+...
        (((1+beta).*(1-(((1-t_d2).^N(m)).*((1-t_d1).^M(m)))))./(2.*(1-(((1-t_d2).^N(m)).*((1-t_d1).^M(m)))+beta)));
    disp('Equilibrium points');
    disp([t_d1 t_d2]);
    disp('Corresponding age and throughput');
    disp([age_1 age_2]);
end
