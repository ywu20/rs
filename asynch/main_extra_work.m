%% 程序描述：引入禁区机制的改进：为了防止陷入局部极大值点，通过判断每一步best周围是否有超过一定数量的nest，若有，将这部分区域设为禁区，
%% 禁止以后飞进去。
%% 预期结果： 找到所有局部极值点


clc;
clear;
global c;
%% 实验的参数设计

Max_condition=3;%考虑Max_condition种不同参数情况，将来准备做两种，第一种是普通cs算法，第二种是引入禁区的改进
L=6;%禁区正方形边长
num_source = 4;
ratio=0.5;%有多少比例要加班
n_tolerance=5;%超过n_tolerance个nest在best附近则判断为禁区
speed=3600;%每个传感器运动速度m/min
N=100;%每种情况运行N次
time=zeros(N,Max_condition);%用来存储Max_condition种不同情况下N次重复实验的时间
successtime0 = 0;
successtime1 = 0;
successtime2=0;
% success_time=zeros(1,Max_condition);%用来记录每一种情况在运行N次中有几次成功找到
%% 算法的参数设计
n=25; %传感器的个数
nd=2; %每个鸟巢的描述维度，具体问题中是每个传感器的二维坐标x，y
pa=0.25; %宿主鸟发现布谷鸟蛋的概率
%% 气体浓度的参数设计
k=0.03;v=0.03;
q0=(2*pi*k)^0.5;x0=10;y0=10; %目标源
q1=0.8*q0;      x1=23;y1=45; %干扰源1
q2=0.7*q0;      x2=45;y2=10; %干扰源2
q3=0.6*q0;      x3=40;y3=45; %干扰源3
q4=0.5*q0;      x4=30;y4=30; %干扰源3
q5=0.8*q0;      x5=13;y5=35; %干扰源1
q6=0.7*q0;      x6=25;y6=16; %干扰源2
q7=0.6*q0;      x7=14;y7=25; %干扰源3
c=gas_concentration_generate(k,v,q0,x0,y0,q1,x1,y1,q2,x2,y2,q3,x3,y3,q4,x4,y4,q5,x5,y5,q6,x6,y6,q7,x7,y7);%多点源产生的浓度空间
%% 算法终止条件
Tmax=500; %最大运行步数
record_max_concentration=zeros(Tmax,N);%用来记录每次运行，每次迭代的最佳浓度，最后根据他画图就ok，如果要展现不同的pa，n等情况，可以把变量扩充为三维。
%% 搜索空间界定
M=50;%搜索方向上限
Lb=1*ones(1,nd);  %x,y方向的最小值
Ub=M*ones(1,nd);  %x,y方向的最大值

%% 传感器初始位置设计
% for i=1:n
%     initial_nest(i,:)=4/5*Ub+Lb+(1/5*Ub-Lb).*rand(size(Lb)); %注：这里通过初始位置设计传感器初始位置不好，离源较远的情况
% end

%% 产生各个鸟巢（传感器）随机的初始位置（Case 1）
for i=1:n
initial_nest(i,:)=Lb+(Ub-Lb).*rand(size(Lb)); %注：这里也可以通过初始位置设计传感器初始位置不好，离源很远的情况
 end

%% 实验开始
for m=1:Max_condition %不同情况（普通cs和改进cs）
    for j=1:N         %每种情况重复进行N次实验
        record_forbidden_center=[];%记录禁区中心坐标\
        forbiddenUb=[];
forbiddenLb=[];
forbidden_center=[];
nest=initial_nest;
   %每次试验传感器初始都一样（控制变量法）
        actual_nest=nest; %记录传感器的真实位置，用于之后计算每次移动的距离distance
        %% 得到当前的最佳传感器
        fitness=10^10*ones(n,1);
        [fmin,bestnest,nest,fitness]=get_best_nest(nest,nest,fitness); %先得到初始时的bestnest等
        
        T=0;
        %% 开始循环迭代
        
        while T<Tmax%这里是采用运行次数到达上界就停止
            T=T+1;
           %% 第一步更新：levy飞行下达指令（最好的sensor不动！！！！）
             new_nest=get_cuckoos(nest,bestnest,Lb,Ub,record_forbidden_center,L);
           
           %% 计算第一步更新的实际上每个sensor的移动距离=levy 飞行指令位置-当前位置
            distance_first=sqrt(sum((new_nest-actual_nest).^2,2));%计算第一步的距离
           %% 进行评估
            [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness);
            actual_nest=new_nest;%更新实际位置
            if m==3  %也可以加上&T<Tmax/2
                new_nest=extra_work(new_nest,best,Lb,Ub,distance_first,speed,ratio);%ratio比例的先到达指令位置的sensor再进行一次限制距离的levy飞行
            else
            %% 第二步更新：pa概率发现并移动
            new_nest=empty_nests(nest,Lb,Ub,pa);%产生第二步移动指令
            end
            % 评估移动
            [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness);

            %% 计算第二步距离以及第一二步总时间
            distance_second=sqrt(sum((new_nest-actual_nest).^2,2));%计算第二步的距离
            actual_nest=new_nest;%更新实际位置
            %sensor是同步的，计算当前一步最长距离，求出该次时长
            distance=distance_first+distance_second;%距离加和得到当前iteration的距离
           if m==1 
               if max(record_max_concentration(:,j)) <= 4.5
            time(j,m)=time(j,m)+max(distance)/speed;%利用最大距离/速度得到截止到当前iteration的总时间
               else
                   time(j,m) = time(j,m);
               end
           end
           %% 在这里引入改进，禁区机制
           if m~=1
            if T>10 %前T/4先不考虑禁区机制
                if fmin~=0
              cuckoo_index=cuckoo_division(nest,best,L);%利用cuckoo_division函数标记在当前best周围的nest
                  if sum(cuckoo_index)>=n_tolerance %如果发现best周围nest超过n_tolerance，则设立禁区
                   local_best=best;
                      forbidden_center = local_best;
                     fUb(1,1) = forbidden_center(1,1)+L/2;
                         fUb(1,2) = forbidden_center(1,2)+L/2;
                         fLb(1,1) = forbidden_center(1,1)-L/2;
                         fLb(1,2) = forbidden_center(1,2)-L/2;
                         a=0;
                         if sum(fLb<Lb)>0 || sum(fUb>Ub)>0
                             continue;
                         end
                        if size(record_forbidden_center,1)~=0
                            record_forbidden_center;
                   for k=1:size(record_forbidden_center,1)
                       overlap = [(forbiddenLb(k,1)<=fLb(1,1) & fLb(1,1)<=forbiddenUb(k,1)) (forbiddenLb(k,2)<=fLb(1,2) & fLb(1,2)<=forbiddenUb(k,2)) (forbiddenLb(k,1)<=fUb(1,1) & fUb(1,1)<=forbiddenUb(k,1)) (forbiddenLb(k,2)<=fUb(1,2) & fUb(1,2)<=forbiddenUb(k,2))];
                       if sum(overlap) >=2
                      forbiddenLb(k,1)=min(forbiddenLb(k,1),fLb(1,1));
                      forbiddenLb(k,2)=min(forbiddenLb(k,2),fLb(1,2));
                      forbiddenUb(k,1)=max(forbiddenUb(k,1),fUb(1,1));
                      forbiddenUb(k,2)=max(forbiddenUb(k,2),fUb(1,2));
                       a=1;
                       end
                   end
                         end  
                   if a==0
                           record_forbidden_center(end+1,:)=local_best;%记录禁区中心坐标
                        forbiddenUb(end+1,:) = fUb;
                         forbiddenLb(end+1,:) = fLb;
                end%记录禁区中心坐标
     end
 end
%                 break;
            
                    
                
              %% 把禁区中的nest踢出
                for num=1:n %每一个nest遍历
                    for k=1:size(record_forbidden_center) %每一个禁区遍历
                     %% 若第num个nest在第k个禁区内，踢出禁区，踢出范围在以禁区中心为中心，边长为5L的正方形（不包括禁区L正方形）范围内
                        if forbiddenLb(k,1)<nest(num,1) & nest(num,1)<forbiddenUb(k,1) & forbiddenLb(k,2)<nest(num,2) & nest(num,2)<forbiddenUb(k,2)
%                             if nest(num,1)<record_forbidden_center(k,1) %若nest在禁区左半边
%                                 nest(num,1)=record_forbidden_center(k,1)-L/1-2*L*rand(1);%则把nest往左边踢出
%                             else
%                                 nest(num,1)=record_forbidden_center(k,1)+L/2+2*L*rand(1);%则把nest往右踢出
%                             end
nest(num,1)=record_forbidden_center(k,2)-5/2*L+5*L*rand(1);%上下方向也要变化，否则levy处陷入死循环
                            nest(num,2)=record_forbidden_center(k,2)-5/2*L+5*L*rand(1);%上下方向也要变化，否则levy处陷入死循环
                        end
                    end
                    nest(num,:)=simplebounds(nest(num,:),Lb,Ub);%要保证踢出后nest在搜索范围内
                    fitness(num)=fobj(nest(num,:));%更新fitness原来禁区里fitness再好也不要了
                end
                [fnew,best,nest,fitness]=get_best_nest(nest,nest,fitness);%驱逐出去后比较一下找best
            end
if size (record_forbidden_center,1) ~= num_source 
    time(j,m)=time(j,m)+max(distance)/speed;
else
    time(j,m) = time(j,m);
end
           end
            %% 找到截止目前的最佳位置bestnest和fmin函数
            if fnew<fmin,
                fmin=fnew;
                bestnest=best;
            end
            record_max_concentration(T,j)=-fmin;%用来存储每一次实验每一次iteration下的最大浓度值，方便之后画图

%             %% 循环截止条件：有1个sensor到目标源附近1m范围内就停止循环
%             if ((bestnest-[x0,y0]).^2)<1
%                 success_time(1,m)=success_time(1,m)+1;%成功找到，成功次数+1
%                 break;
%             end
%             figure(1);
%             hold off;
%             for i=1:n
%                 axis([0 M 0 M]);
%                 plot(nest(i,1),actual_nest(i,2),'*');%%n个备选解（nest）用不同颜色的*表示
%                 plot(best(1,1),best(1,2),'.','markersize',20);%%最优解（bestnest）用加粗的实心点表示
%                 hold on;
%             end
%             x=1:M;
%             y=1:M;
%             global c;
%             contour(x,y,c);
%             for i=1:size(record_forbidden_center,1)
%                 hold on
%                 center=record_forbidden_center(i,:);%中心
%                 length=sqrt(2)/2*L; %边长
%                 rotate=0*pi; %旋转角度
%                 aa=[-1 -1 1 1 -1];
%                 bb=[-1 1 1 -1 -1];
%                 cc=complex(aa,bb)/sqrt(2);
%                 cc=cc*exp(j*rotate)*length+complex(center(1),center(2));
%                 line(real(cc),imag(cc),'LineWidth',1);
%             end
%             title(num2str(T));
% %             pause(1);
%         end
%          if m==2
%         figure(4)
% t=1:Tmax;
% plot (t,size(record_forbidden_center,1))
          end% iteration结束
          if max(record_max_concentration(:,j))>= 4 & m==1
%               time (j,m) = 1000;
 successtime0 = successtime0 + 1;
          elseif size(record_forbidden_center,1) == num_source & m==2
%               time(j,m) = 1000;
successtime1 = successtime1 + 1;
  elseif size(record_forbidden_center,1) == num_source & m==3
%               time(j,m) = 1000;
successtime2 = successtime2 + 1;
          end
     end     % 重复实验结束
end   % 不同情况结束 
% 
% bestnest %在command面板上展示出来
% % figure (2)
% % boxplot (record_max_concentration)
% for i=1:size(record_forbidden_center,1)
%     hold on;
%     center=record_forbidden_center(i,:);%中心
%     length=sqrt(2)/2*L; %边长
%     rotate=0*pi; %旋转角度
%     aa=[-1 -1 1 1 -1];
%     bb=[-1 1 1 -1 -1];
%     cc=complex(aa,bb)/sqrt(2);
%     cc=cc*exp(j*rotate)*length+complex(center(1),center(2));
%     line(real(cc),imag(cc),'LineWidth',1);
% end
% %     axis equal;
record_forbidden_center
% 三个评价指标：
% success_time %成功找到的次数对比
% sum(time(:,1)>time(:,2)) %N次实验中，加班改进后用时少的次数
% average_time=mean(time) %N此实验平均用时对比 
figure(4) 
successstat = [successtime0 successtime1 successtime2];
X=1:1:3;
bar(X,successstat);
title('Compare the number of times the requirement is met in 100 generations');
%% 画图

figure (3)
%average_time=mean(time);%求N次重复实验的平均时间
%plot(1:Max_condition,average_time,'.','MarkerSize',30)
boxplot (time, 'Labels', {'conventional','improved', 'extra-improved'});
title('Compare the time it takes for the robots to meet the requirements');
