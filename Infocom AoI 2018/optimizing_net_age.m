function age = optimizing_net_age(beta,t_d,t_w,min,max,M,N)
a = 0;
b = 1;
age = ((1-(((1-t_d)^M)*((1-t_w)^N))+beta)/(t_d*((1-t_d)^(M-1))*((1-t_w)^N)))+(beta/2)+...
    (((1+beta)*(1-(((1-t_d)^M)*((1-t_w)^N))))/(2*(1-(((1-t_d)^M)*((1-t_w)^N))+beta)));
age = (((b-a)*(age-min))/(max - min)) + a;
cost = 0;
age = age + cost;