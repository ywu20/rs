function z=fobj(u)
x=u(1,1);
y=u(1,2);
global c
c_measure=c(floor(y),floor(x));

z=-c_measure;
