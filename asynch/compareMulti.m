%asynchronous cuckoo search for finding one absolute maximum
%function compare(initial_nest)
global c;
%% 实验的参数设计
speed=3600;%每个传感器运动速度m/min
N=100;%每种情况运行N次
time=zeros(N);%用来存储Max_condition种不同情况下N次重复实验的时间
successtime0 = 0; %记录成功次数
successtime1=0;
L=6;%禁区正方形边长
num_source = 4;
n_tolerance=5;%超过n_tolerance个nest在best附近则判断为禁区
record_forbidden_center=[];
%% 算法的参数设计
n=50; %传感器的个数
nd=2; %每个鸟巢的描述维度，具体问题中是每个传感器的二维坐标x，y
pa=0.25; %宿主鸟发现布谷鸟蛋的概率
cuckooInfo = zeros(n,8);
Tmax=3000;
%% 气体浓度的参数设计
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

%% 搜索空间界定
M=50;%搜索方向上限
Lb=1*ones(1,nd);  %x,y方向的最小值
Ub=M*ones(1,nd);  %x,y方向的最大值

%% far from source distribution
% for i=1:n
%     initial_nest(i,:)=4/5*Ub+Lb+(1/5*Ub-Lb).*rand(size(Lb)); %注：这里通过初始位置设计传感器初始位置不好，离源较远的情况
% end

%%uniform random distribution
% 
for i=1:n
initial_nest(i,:)=Lb+(Ub-Lb).*rand(size(Lb)); %注：这里也可以通过初始位置设计传感器初始位置不好，离源很远的情况 
end
successstat = zeros(1,2);
Max_condition=2;
time1 = zeros(1,N);
time2 = zeros(1,N);

 for m=1:Max_condition %不同情况（普通cs和改进cs）
for j=1:N
   if m==1
       [successtime0 time1(1,j)] = main_2(initial_nest);
       if(successtime0 == 1)
           successstat(1,1) = successstat(1,1)+1;
        end 
   else
        [successtime1 time2(1,j)]=multi_source(initial_nest);
   if(successtime1==1)
       successstat(1,2) = successstat(1,2) +1;
   end
        
   end
end
 end

figure(4) 
  X = categorical({'classical','multiSource'});
  successstat
 bar(X,successstat);
 title('Accuracy of Two Algorithms');
set(gca,'FontSize',14,'Fontname', 'Times New Roman');
% %  % %% 画图
% % % 
% % 
figure (3)
TIME=zeros(N,Max_condition);
TIME(:,1)=time1;%求N次重复实验的平均时间
TIME(:,2) = time2;
% %plot(1:Max_condition,average_time,'.','MarkerSize',30)
boxplot (TIME, 'Labels', {'classical','multiSource'});
set(gca,'FontSize',14,'Fontname', 'Times New Roman');
median(TIME(:,1))
median(TIME(:,2))
title('Searching Time');
%end