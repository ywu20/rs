function c=gas_concentration_generate(k,v,q0,x0,y0,q1,x1,y1,q2,x2,y2,q3,x3,y3,q4,x4,y4,q5,x5,y5,q6,x6,y6,q7,x7,y7)

x=1:1:60;
y=1:1:60;
[x y]=meshgrid(x,y);
d0=sqrt((x-x0).^2+(y-y0).^2);
d0(x0,y0)=1;
d1=sqrt((x-x1).^2+(y-y1).^2);
d1(y1,x1)=1;
d2=sqrt((x-x2).^2+(y-y2).^2);
d2(y2,x2)=1;
d3=sqrt((x-x3).^2+(y-y3).^2);
d3(y3,x3)=1;
d4=sqrt((x-x4).^2+(y-y4).^2);
d4(y4,x4)=1;
d5=sqrt((x-x5).^2+(y-y5).^2);
d5(y5,x5)=1;
d6=sqrt((x-x6).^2+(y-y6).^2);
d6(y6,x6)=1;
d7=sqrt((x-x7).^2+(y-y7).^2);
d7(x7,y7)=1;

c0=q0*exp(-v/(2*k).*(d0-(x-x0)))./(2*pi*k.*d0);
c0(y0,x0)=4.6;
c1=q1*exp(-v/(2*k).*(d1-(x-x1)))./(2*pi*k.*d1);
c1(y1,x1)=3.6;
c2=q2*exp(-v/(2*k).*(d2-(x-x2)))./(2*pi*k.*d2);
c2(y2,x2)=3.2;
c3=q3*exp(-v/(2*k).*(d3-(x-x3)))./(2*pi*k.*d3);
c3(y3,x3)=2.6;
c4=q4*exp(-v/(2*k).*(d4-(x-x4)))./(2*pi*k.*d4);
c4(y4,x4)=2.6;
c5=q5*exp(-v/(2*k).*(d5-(x-x5)))./(2*pi*k.*d5);
c5(y5,x5)=2.7;
c6=q6*exp(-v/(2*k).*(d6-(x-x6)))./(2*pi*k.*d6);
c6(y6,x6)=3.2;
c7=q7*exp(-v/(2*k).*(d7-(x-x7)))./(2*pi*k.*d7);
c7(y7,x7)=2.6;
%
c=c0+c1+c2+c3;%+c4+c5+c6+c7;

% figure (4)
% colormap(prism)
% surf(c)
% figure(5);
% contour(x,y,c);
% axis([0 50 0 50]);
