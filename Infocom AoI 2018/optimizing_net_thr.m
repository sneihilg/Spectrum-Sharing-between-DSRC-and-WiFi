function thr = optimizing_net_thr(beta,t_w,t_d,M,N)
thr = -(N*t_w*((1-t_w)^(N-1))*((1-t_d)^M)*(1+beta))/(1-(((1-t_d)^M)*((1-t_w)^N))+beta);
cost = 0;
thr = thr + cost;