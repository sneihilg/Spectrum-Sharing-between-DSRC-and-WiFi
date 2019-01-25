function thr = stackelberg_optimize_net_thr(beta,t_w,w_col,w_idle,M,N,options)
% find optimal t_d
t_d_o = 1e-2;
lb = 1e-2;
ub = 0.99;
t_d_temp = 1e-2:1e-2:0.99;
t_w_temp = t_w;
age_temp = ((((1-(((1-t_d_temp).^M).*((1-t_w_temp).^N))+beta)./(t_d_temp.*((1-t_d_temp).^(M-1)).*((1-t_w_temp).^N)))+(beta/2)+...
        (((1+beta)*(1-(((1-t_d_temp).^M).*((1-t_w_temp).^N))))./(2*(1-(((1-t_d_temp).^M).*((1-t_w_temp).^N))+beta)))));
minimum = min(age_temp);
maximum = max(age_temp);
t_d = fmincon(@(t_d)optimizing_net_age_cost(beta,t_d,t_w,w_col,w_idle,minimum,maximum,M,N),t_d_o,[],[],[],[],lb,ub,[],options);

thr = -(N*t_w*((1-t_w)^(N-1))*((1-t_d)^M)*(1+beta))/(1-(((1-t_d)^M)*((1-t_w)^N))+beta)+...
    w_idle*(((1-t_d)^M)*((1-t_w)^N)) + w_col*(1-(((1-t_d)^M)*((1-t_w)^N))...
    - (M*t_d*((1-t_d)^(M-1))*((1-t_w)^N)) - (N*t_w*((1-t_w)^(N-1))*((1-t_d)^M)));