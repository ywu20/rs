global c;
Max_condition =4;
N=10;
n=30;
nd=2;
M=50; %������������
Lb=1*ones(1,nd);   %x,y�������Сֵ
Ub=M*ones(1,nd);
k=0.03;v=0.03;
q0=(2*pi*k)^0.5;x0=10;y0=10; %Ŀ��Դ
q1=0.8*q0;      x1=23;y1=45; %����Դ1
q2=0.7*q0;      x2=45;y2=10; %����Դ2
q3=0.6*q0;      x3=40;y3=45; %����Դ3
q4=0.5*q0;      x4=30;y4=30; %����Դ4
q5=0.8*q0;      x5=13;y5=35; %����Դ5
q6=0.7*q0;      x6=25;y6=16; %����Դ6
q7=0.6*q0;      x7=14;y7=25; %����Դ7
c=gas_concentration_generate(k,v,q0,x0,y0,q1,x1,y1,q2,x2,y2,q3,x3,y3,q4,x4,y4,q5,x5,y5,q6,x6,y6,q7,x7,y7);%���Դ������Ũ�ȿռ�
 successtime1=0;
 successtime2 = 0;successtime0=0;
successtime3=0;

successstat = zeros(1,2);
Max_condition=4;
time1 = zeros(1,N);
time2 = zeros(1,N);
time3 = zeros(1,N);
time4 = zeros(1,N);
doesntwork = 0; dw =0;
 for m=1:Max_condition %��ͬ�������ͨcs�͸Ľ�cs��
for j=1:N
    for i=1:n
    initial_nest(i,:)=Lb+(Ub-Lb).*rand(size(Lb)); %ע������Ҳ����ͨ����ʼλ����ƴ�������ʼλ�ò��ã���Դ��Զ����� 
    end
   if m==1
       %successtime0=0;
       %while successtime0 == 0
           [successtime time1(1,j)] = asynchronous_cuckoo_search(initial_nest);
           dw=dw+1;
           successtime0 = successtime0 +successtime;
       %end
   elseif m==2
      
      % while successtime1 == 0
           [successtime time2(1,j)] = main_2(initial_nest);
      % end
      successtime1 = successtime1 +successtime;
   elseif m== 3
       
       %while successtime2 == 0
           doesntwork=doesntwork+1;
           [successtime time3(1,j)] = asynch_pso(initial_nest);
           successtime2 = successtime2 +successtime;
      % end
       
   else 
       
       %while successtime3 == 0
       [successtime time4(1,j)] = pso(initial_nest);
       successtime3 = successtime3+successtime;
       %end
   end
end
 end
successtime0
successtime1
successtime2
successtime3
r = randi([-10 10],1,1000)
%%
% 
%   for x = 1:10
%       disp(x)
%   end
% 
 figure (3)
TIME=zeros(N,2);
TIME(:,1)=time1;%��N���ظ�ʵ���ƽ��ʱ��
TIME(:,2) = time2;
%TIME(:,3)=time3;%��N���ظ�ʵ���ƽ��ʱ��
%TIME(:,4) = time4;

% %plot(1:Max_condition,average_time,'.','MarkerSize',30)
boxplot (TIME, 'Labels', {'asynchronous csa','classical csa'});
set(gca,'FontSize',14,'Fontname', 'Times New Roman');
median(TIME(:,1))
median(TIME(:,2))
%median(TIME(:,3))
%median(TIME(:,4))
title('Searching Time');

figure (4)
TIME4=zeros(N,2);
TIME4(:,1)=time1;%��N���ظ�ʵ���ƽ��ʱ��
TIME4(:,2) = time3;
%TIME(:,3)=time3;%��N���ظ�ʵ���ƽ��ʱ��
%TIME(:,4) = time4;

% %plot(1:Max_condition,average_time,'.','MarkerSize',30)
boxplot (TIME4, 'Labels', {'asynchronous csa', 'asynchronous pso'});
set(gca,'FontSize',14,'Fontname', 'Times New Roman');
median(TIME4(:,1))
median(TIME4(:,2))
%median(TIME(:,3))
%median(TIME(:,4))
title('Searching Time');