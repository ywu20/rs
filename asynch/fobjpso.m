function z=fobjpso(u,c)
%% d-dimensional sphere function sum_j=1^d (u_j-1)^2.
% with a minimum at (1,1, ...., 1);
%%u�����зֱ�Ϊʵ�ʴ�������λ�������x,y
x=u(1,1)
y=u(1,2)

c_measure=c(floor(y),floor(x));

z=-c_measure;%Ϊ��������Сֵ������ȡ�෴��
