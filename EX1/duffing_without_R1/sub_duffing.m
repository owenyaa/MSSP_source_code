function yprime=sub_duffing(t,y)
global A epsilon
yprime=A*y+[0;-epsilon].*[0;y(1)^3];