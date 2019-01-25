function age = stackelberg_optimize_net_age(beta,t_d,w_col,w_idle,min,max,M,N,options)
a = 0;
b = 1;
% find WiFi's best response to DSRC node
t_w_o = 0.99;
lb = 1e-2;
ub = 0.99;
t_w = fmincon(@(t_w)optimizing_net_thr_cost(beta,t_w,t_d,w_col,w_idle,M,N),t_w_o,[],[],[],[],lb,ub,[],options);

age = ((1-(((1-t_d)^M)*((1-t_w)^N))+beta)/(t_d*((1-t_d)^(M-1))*((1-t_w)^N)))+(beta/2)+...
    (((1+beta)*(1-(((1-t_d)^M)*((1-t_w)^N))))/(2*(1-(((1-t_d)^M)*((1-t_w)^N))+beta)));
age = (((b-a)*(age-min))/(max - min)) + a;
cost = w_idle*(((1-t_d)^M)*((1-t_w)^N)) + w_col*(1-(((1-t_d)^M)*((1-t_w)^N))...
    -(M*t_d*((1-t_d)^(M-1))*((1-t_w)^N)) - (N*t_w*((1-t_w)^(N-1))*((1-t_d)^M)));
age = age + cost;