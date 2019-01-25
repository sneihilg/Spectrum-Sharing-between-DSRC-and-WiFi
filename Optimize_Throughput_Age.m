% code to generate the value of optimal tau(s) that maximizes (a) the
% throughput of nodes in 2 WiFi networks (b) that minimizes the age of
% nodes in 2 DSRC networks

clc;
clear all;
close all;

beta = 0.001;     
l_idle = beta;          %length of an idle slot
l_col = 1+beta;         %length of collision slot

M = [1 2 5];
N = [1 2 5];

n = M+N;

options = optimset('Algorithm','interior-point','TolX',1e-14,...
     'TolFun',1e-12,'TolCon',1e-14,'MaxFunEval',1e8,'MaxIter',1e8,'Display','Iter','InitBarrierParam',1e-18);

optimal_tau_w = [];
optimal_thr = [];
optimal_tau_d = [];
optimal_age = [];
for m = 1:numel(M)
    t_w_o = 1e-2;
    lb = 1e-2;
    ub = 0.99;
    [t_w,fval_thr,exitflag] = fmincon(@(t_w)optimize_throughput(l_col,l_idle,t_w,n(m)),t_w_o,[],[],[],[],lb,ub,[],options);  
    optimal_tau_w = [optimal_tau_w t_w];
    optimal_thr = [optimal_thr -fval_thr/n(m)];
    
    t_d_o = 1e-2;
    lb = 1e-2;
    ub = 0.99;
    [t_d,fval_age,exitflag] = fmincon(@(t_d)optimize_age(beta,t_d,n(m)),t_d_o,[],[],[],[],lb,ub,[],options);
    optimal_tau_d = [optimal_tau_d t_d];
    optimal_age = [optimal_age fval_age/n(m)];
end
disp('Optimal access probability (WiFi nodes)');
disp(optimal_tau_w);
disp('Optimal throughput single var');
disp(optimal_thr);
disp('Optimal access probability (DSRC nodes)');
disp(optimal_tau_d);
disp('Optimal age');
disp(optimal_age);