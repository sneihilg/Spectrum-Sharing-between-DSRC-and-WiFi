function thr = optimize_throughput(l_col,l_idle,t_w,n)
thr = -(n*t_w*((1-t_w)^(n-1))*l_col)/(1-((1-t_w)^n) + l_idle);