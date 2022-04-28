
co_z = [0 0 0.08 0.08 0.04 0 0];

t = linspace(0,1,101);
y = polyval_bz(co_z,t);
dy = polyval_bz_d(co_z,t);

yyaxis left
plot(t,y)

yyaxis right
plot(t,dy)




