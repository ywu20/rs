global c;
Max_condition =4;
N=1;
n=30;
nd=2;
M=50; %搜索方向上限
Lb=1*ones(1,nd);   %x,y方向的最小值
Ub=M*ones(1,nd);
k=0.03;v=0.03;
q0=(2*pi*k)^0.5;x0=10;y0=10; %目标源
q1=0.8*q0;      x1=23;y1=45; %干扰源1
q2=0.7*q0;      x2=45;y2=10; %干扰源2
q3=0.6*q0;      x3=40;y3=45; %干扰源3
q4=0.5*q0;      x4=30;y4=30; %干扰源4
q5=0.8*q0;      x5=13;y5=35; %干扰源5
q6=0.7*q0;      x6=25;y6=16; %干扰源6
q7=0.6*q0;      x7=14;y7=25; %干扰源7
c=gas_concentration_generate(k,v,q0,x0,y0,q1,x1,y1,q2,x2,y2,q3,x3,y3,q4,x4,y4,q5,x5,y5,q6,x6,y6,q7,x7,y7);%多点源产生的浓度空间

[successtime time] = asynch_pso(initial_nest);

successtime
time

