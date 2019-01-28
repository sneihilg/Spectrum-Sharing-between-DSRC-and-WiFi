function age = optimizing_net_age_cost(beta,t_d,t_w,w_col,w_idle,min,max,M,N)
a = 0;
b = 1;
age = ((1-(((1-t_d)^M)*((1-t_w)^N))+beta)/(t_d*((1-t_d)^(M-1))*((1-t_w)^N)))+(beta/2)+...
    (((1+beta)*(1-(((1-t_d)^M)*((1-t_w)^N))))/(2*(1-(((1-t_d)^M)*((1-t_w)^N))+beta)));
age = (((b-a)*(age-min))/(max - min)) + a;
cost = w_idle*(((1-t_d)^M)*((1-t_w)^N)) + w_col*(1-(((1-t_d)^M)*((1-t_w)^N))...
    -(M*t_d*((1-t_d)^(M-1))*((1-t_w)^N)) - (N*t_w*((1-t_w)^(N-1))*((1-t_d)^M)));
age = age + cost;