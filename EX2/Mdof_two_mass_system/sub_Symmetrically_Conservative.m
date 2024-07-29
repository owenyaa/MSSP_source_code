function yprime=sub_Symmetrically_Conservative(t,y)
global A k22
yprime=A*y-k22*[0;0;(y(1)-y(2))^3;(y(2)-y(1))^3];