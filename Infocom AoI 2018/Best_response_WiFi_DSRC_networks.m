% plotting best response curve for age and throughput nodes without cost
clc;
clear all;
close all;

beta = 0.001;     
l_idle = beta;          %length of an idle slot
l_col = 1+beta;         %length of collision slot
M = [5 5 5];                  % number of DSRC nodes
N = [1 2 5];                  % number of WiFi nodes
tau_w = 1e-2:1e-2:0.99;
tau_d = 1e-2:1e-2:0.99;

if M(1)<5
    options = optimset('Algorithm','interior-point','TolX',1e-14,...
         'TolFun',1e-12,'TolCon',1e-14,'MaxFunEval',1e6,'MaxIter',1e6,'Display','Iter','InitBarrierParam',1e-12);
else
    options = optimset('Algorithm','interior-point','TolX',1e-12,...
     'TolFun',1e-10,'TolCon',1e-12,'MaxFunEval',1e6,'MaxIter',1e6,'Display','Iter','InitBarrierParam',1e-2);
end

optimal_tau_w = [];
optimal_tau_d = [];
NE = [];
for m = 1:numel(M)
    optim_tau_w = [];
    optim_thr = [];
    for i = 1:numel(tau_d)
        t_w_o = 1e-2;
        lb = 1e-2;
        ub = 0.99;
        [t_w,fval,exitflag] = fmincon(@(t_w)optimizing_net_thr(beta,t_w,tau_d(i),M(m),N(m)),t_w_o,[],[],[],[],lb,ub,[],options);
        optim_tau_w = [optim_tau_w t_w];
        optim_thr = [optim_thr -fval];
    end
    optimal_tau_w = [optimal_tau_w ; optim_tau_w];
    
    optim_tau_d = [];
    optim_age = [];
    for i = 1:numel(tau_w)
        t_d_o = 1e-2;
        lb = 1e-2;
        ub = 0.99;
        t_d_temp = 1e-2:1e-2:0.99;
        t_w_temp = tau_w(i);
        age_temp = ((((1-(((1-t_d_temp).^M(m)).*((1-t_w_temp).^N(m)))+beta)./(t_d_temp.*((1-t_d_temp).^(M(m)-1)).*((1-t_w_temp).^N(m))))+(beta/2)+...
                (((1+beta)*(1-(((1-t_d_temp).^M(m)).*((1-t_w_temp).^N(m)))))./(2*(1-(((1-t_d_temp).^M(m)).*((1-t_w_temp).^N(m)))+beta)))));
        minimum = min(age_temp);
        maximum = max(age_temp);
        [t_d,fval,exitflag] = fmincon(@(t_d)optimizing_net_age(beta,t_d,tau_w(i),minimum,maximum,M(m),N(m)),t_d_o,[],[],[],[],lb,ub,[],options);
        optim_tau_d = [optim_tau_d t_d];
        optim_age = [optim_age fval];
    end
    optimal_tau_d = [optimal_tau_d ; optim_tau_d];

    % find the point of intersection on the plot
    [xout,yout] = intersections(optim_tau_d,tau_w,tau_d,optim_tau_w,1);
    NE{m} = [xout yout];
    
    figure;
    plot(optimal_tau_d(m,:),tau_w,'-r','LineWidth',8,'MarkerSize',18);
    hold on;
    plot(tau_d,optimal_tau_w(m,:),':k','LineWidth',8,'MarkerSize',18);
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
    fig_name = sprintf('M_%d_N_%d.pdf',M(m),N(m));
    xlabel('\tau_D','fontweight','bold','fontsize',42);
    ylabel('\tau_W','fontweight','bold','fontsize',42);
    if M(m)<5
        legend('DSRC','WiFi','Location','northwest');
    else
        legend('DSRC','WiFi','Location','northeast');
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
    t_d = NE{m}(:,1);
    t_w = NE{m}(:,2);
    thr = (t_w.*((1-t_w).^(N(m)-1)).*((1-t_d).^M(m)).*(1+beta))./(1-(((1-t_d).^M(m)).*((1-t_w).^N(m)))+beta);
    age = ((1-(((1-t_d).^M(m)).*((1-t_w).^N(m)))+beta)./(t_d.*((1-t_d).^(M(m)-1)).*((1-t_w).^N(m))))+(beta/2)+...
        (((1+beta).*(1-(((1-t_d).^M(m)).*((1-t_w).^N(m)))))./(2.*(1-(((1-t_d).^M(m)).*((1-t_w).^N(m)))+beta)));
    disp('Equilibrium points');
    disp([t_d t_w]);
    disp('Corresponding age and throughput');
    disp([age thr]);
end
