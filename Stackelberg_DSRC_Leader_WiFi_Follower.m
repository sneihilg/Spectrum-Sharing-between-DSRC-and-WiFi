% Stackelberg game for two network game
clc;
close all;
clear all;

beta = 0.001;     
l_idle = beta;          %length of an idle slot
l_col = 1+beta;         %length of collision slot
M = 5;                  % number of DSRC nodes
N = 5;                  % number of WiFi nodes
w_col = 1+beta;
w_idle = beta;

options = optimset('Algorithm','interior-point','TolX',1e-14,...
     'TolFun',1e-12,'TolCon',1e-14,'MaxFunEval',1e6,'MaxIter',1e6,'Display','Iter','InitBarrierParam',1e-18);

optimal_tau_w = [];
optimal_tau_d = [];
SE = [];
for m = 1:numel(M) 
    % find the value of \tau_d for the DSRC node
    t_d_temp = 1e-2:1e-2:0.99;
    t_w_temp = [];
    for i=1:numel(t_d_temp)
        t_w_o = 1e-2;
        lb = 1e-2;
        ub = 0.99;
        [t_w,fval,exitflag] = fmincon(@(t_w)optimizing_net_thr_cost(beta,t_w,t_d_temp(i),w_col,w_idle,M(m),N(m)),t_w_o,[],[],[],[],lb,ub,[],options);
        t_w_temp = [t_w_temp t_w];
    end

    t_d_o = 1e-2;
    age_temp = (((1-((1-t_d_temp).*(1-t_w_temp))+beta)./(t_d_temp.*(1-t_w_temp)))+(beta/2)+...
    (((1+beta)*(1-((1-t_d_temp).*(1-t_w_temp))))./(2*(1-((1-t_d_temp).*(1-t_w_temp))+beta))));
    minimum = min(age_temp);
    maximum = max(age_temp);
    lb = 1e-2;
    ub = 0.99;
    [t_d,fval,exitflag] = fmincon(@(t_d)stackelberg_optimize_net_age(beta,t_d,w_col,w_idle,minimum,maximum,M(m),N(m),options),t_d_o,[],[],[],[],lb,ub,[],options); 

    % find the best response of WiFi node to the DSRC node's strategy
    t_w_o = 1e-2;
    lb = 1e-2;
    ub = 0.99;
    [t_w,fval,exitflag] = fmincon(@(t_w)optimizing_net_thr_cost(beta,t_w,t_d,w_col,w_idle,M(m),N(m)),t_w_o,[],[],[],[],lb,ub,[],options);
    
    SE = [SE ; t_d t_w];
end

for m = 1:numel(M)
    % final computation
    t_d = SE(m,1);
    t_w = SE(m,2);
    thr = (t_w.*((1-t_w).^(N(m)-1)).*((1-t_d).^M(m)).*(1+beta))./(1-(((1-t_d).^M(m)).*((1-t_w).^N(m)))+beta);
    cost_thr = w_idle.*(((1-t_d).^M(m)).*((1-t_w).^N(m))) + w_col.*(1-(((1-t_d).^M(m)).*((1-t_w).^N(m)))...
        - (M(m).*t_d.*((1-t_d).^(M(m)-1)).*((1-t_w).^N(m))) - (N(m).*t_w.*((1-t_w).^(N(m)-1)).*((1-t_d).^M(m))));
    thr_payoff = thr - cost_thr;
    age = ((1-(((1-t_d).^M(m)).*((1-t_w).^N(m)))+beta)./(t_d.*((1-t_d).^(M(m)-1)).*((1-t_w).^N(m))))+(beta/2)+...
        (((1+beta).*(1-(((1-t_d).^M(m)).*((1-t_w).^N(m)))))./(2.*(1-(((1-t_d).^M(m)).*((1-t_w).^N(m)))+beta)));
    cost_age = w_idle.*(((1-t_d).^M(m)).*((1-t_w).^N(m))) + ...
        w_col.*(1-(((1-t_d).^M(m)).*((1-t_w).^N(m)))-(M(m).*t_d.*((1-t_d).^(M(m)-1)).*((1-t_w).^N(m)))...
        - (N(m).*t_w.*((1-t_w).^(N(m)-1)).*((1-t_d).^M(m))));
    age_payoff = age + cost_age;
    disp('Equilibrium points');
    disp([t_d t_w]);
    disp('Corresponding age and throughput');
    disp([age thr]);
    disp('Corresponding age and throughput payoff');
    disp([age_payoff thr_payoff]);
end