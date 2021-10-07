function z=fobjpso(u,c)
%% d-dimensional sphere function sum_j=1^d (u_j-1)^2.
% with a minimum at (1,1, ...., 1);
%%u的三行分别为实际传感器的位置坐标的x,y
x=u(1,1)
y=u(1,2)

c_measure=c(floor(y),floor(x));

z=-c_measure;%为了搜索最小值，所以取相反数
