global c;

N=1000; % number of iteration being ran
n=30; % number of robots
nd=2; % dimension
M=50; %bound
Lb=1*ones(1,nd);   %lower bound
Ub=M*ones(1,nd); %upper bound
k=0.03;v=0.03; % gas constants
q0=(2*pi*k)^0.5;x0=10;y0=10; 
q1=0.8*q0;      x1=23;y1=45; 
q2=0.7*q0;      x2=45;y2=10; 
q3=0.6*q0;      x3=40;y3=45; 
q4=0.5*q0;      x4=30;y4=30; 
q5=0.8*q0;      x5=13;y5=35; 
q6=0.7*q0;      x6=25;y6=16; 
q7=0.6*q0;      x7=14;y7=25; 
c=gas_concentration_generate(k,v,q0,x0,y0,q1,x1,y1,q2,x2,y2,q3,x3,y3,q4,x4,y4,q5,x5,y5,q6,x6,y6,q7,x7,y7);%多点源产生的浓度空间
 successtime1=0;
 successtime2 = 0;successtime0=0;
successtime3=0;

successstat = zeros(1,2);
Max_condition=3;
time1 = zeros(1,N);
time2 = zeros(1,N);
time3 = zeros(1,N);
time4 = zeros(1,N);
doesntwork = 0; dw =0;

 for j=1:N %不同情况（普通cs和改进cs）
    for i=1:n
    initial_nest(i,:)=Lb+(Ub-Lb).*rand(size(Lb)); %注：这里也可以通过初始位置设计传感器初始位置不好，离源很远的情况 
    end
for m = 1:Max_condition

   if m==1
       %successtime0=0;
       %while successtime0 == 0
           [successtime time1(1,j)] = asynchronous_cuckoo_search(initial_nest);
           successtime0 = successtime0 +successtime;
       %end
   elseif m==2
      
      % while successtime1 == 0
           [successtime time2(1,j)] = main_2(initial_nest);
      % end
      successtime1 = successtime1 +successtime;
   elseif m== 3
       
       %while successtime2 == 0
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
%r = randi([-10 10],1,1000)
%%
% 
%   for x = 1:10
%       disp(x)
%   end
% 
 figure (3)
TIME=zeros(N,2);
TIME(:,1)=time1;%求N次重复实验的平均时间
TIME(:,2) = time2;
%TIME(:,3)=time3;%求N次重复实验的平均时间
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
TIME4(:,1)=time1;%求N次重复实验的平均时间
TIME4(:,2) = time3;
%TIME(:,3)=time3;%求N次重复实验的平均时间
%TIME(:,4) = time4;

% %plot(1:Max_condition,average_time,'.','MarkerSize',30)
boxplot (TIME4, 'Labels', {'asynchronous csa', 'asynchronous pso'});
set(gca,'FontSize',14,'Fontname', 'Times New Roman');
median(TIME4(:,1))
median(TIME4(:,2))
%median(TIME(:,3))
%median(TIME(:,4))
title('Searching Time');