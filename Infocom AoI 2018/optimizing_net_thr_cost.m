function thr = optimizing_net_thr_cost(beta,t_w,t_d,w_col,w_idle,M,N)
thr = -(N*t_w*((1-t_w)^(N-1))*((1-t_d)^M)*(1+beta))/(1-(((1-t_d)^M)*((1-t_w)^N))+beta);
cost = w_idle*(((1-t_d)^M)*((1-t_w)^N)) + w_col*(1-(((1-t_d)^M)*((1-t_w)^N))...
     - (M*t_d*((1-t_d)^(M-1))*((1-t_w)^N)) - (N*t_w*((1-t_w)^(N-1))*((1-t_d)^M)));
thr = thr + cost;