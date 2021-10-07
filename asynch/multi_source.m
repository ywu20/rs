%% 程序描述：引入禁区机制的改进：为了防止陷入局部极大值点，通过判断每一步best周围是否有超过一定数量的nest，若有，将这部分区域设为禁区，
%% 禁止以后飞进去。
%% 预期结果： 找到所有局部极值点


function [success time] = multi_source(initial_nest)
L=6;%禁区正方形边长

n_tolerance=5;%超过n_tolerance个nest在best附近则判断为禁区
speed=3600;%每个传感器运动速度m/min
time=0;%用来存储Max_condition种不同情况下N次重复实验的时间
n=size(initial_nest, 1); %传感器的个数
nd=2; %每个鸟巢的描述维度，具体问题中是每个传感器的二维坐标x，y
pa=0.25; %宿主鸟发现布谷鸟蛋的概率
num_source = 4;
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
%% 搜索空间界定
M=50;%搜索方向上限
Lb=1*ones(1,nd);  %x,y方向的最小值
Ub=M*ones(1,nd);  %x,y方向的最大值
T=0;
        record_forbidden_center=[];%记录禁区中心坐标\
        forbiddenUb=[];
forbiddenLb=[];
forbidden_center=[];
nest=initial_nest;
   %每次试验传感器初始都一样（控制变量法）
        actual_nest=nest; %记录传感器的真实位置，用于之后计算每次移动的距离distance
        %% 得到当前的最佳传感器
        fitness=10^10*ones(n,1);
        [fmin,bestnest,nest,fitness]=get_best_nest2(nest,nest,fitness); %先得到初始时的bestnest等
        %% 开始循环迭代
tStart=tic;
tEnd=0;
        while (tEnd<0.15)%这里是采用运行次数到达上界就停止
     T=T+1;
           %% 第一步更新：levy飞行下达指令（最好的sensor不动！！！！）
            new_nest=get_cuckoos(nest,bestnest,Lb,Ub,record_forbidden_center,L); 
           %% 计算第一步更新的实际上每个sensor的移动距离=levy 飞行指令位置-当前位置
            distance_first=sqrt(sum((new_nest-actual_nest).^2,2));%计算第一步的距离
           %% 进行评估
            [fnew,best,nest,fitness]=get_best_nest2(nest,new_nest,fitness);
            actual_nest=new_nest;%更新实际位置

            %% 第二步更新：pa概率发现并移动
            new_nest=empty_nests(nest,Lb,Ub,pa);%产生第二步移动指令
            % 评估移动
            [fnew,best,nest,fitness]=get_best_nest2(nest,new_nest,fitness);

            %% 计算第二步距离以及第一二步总时间
            distance_second=sqrt(sum((new_nest-actual_nest).^2,2));%计算第二步的距离
            actual_nest=new_nest;%更新实际位置
            %sensor是同步的，计算当前一步最长距离，求出该次时长
            distance=distance_first+distance_second;%距离加和得到当前iteration的距离
            done = 0;
            if(size(record_forbidden_center,1)==num_source)
            for i=1:size(record_forbidden_center,1)
                if sum(abs(bestnest-[x0 y0]))<=5||sum(abs(bestnest-[x1 y1]))<= 5 ||sum(abs(bestnest-[x2 y2]))<= 5 || sum(abs(bestnest-[x3 y3]))<=5 || sum(abs(bestnest-[x4 y4]))<=5||sum(abs(bestnest-[x5 y5]))<= 5 ||sum(abs(bestnest-[x6 y6]))<= 5    
                   done=done+1;
                end
            end
            end
              if done ~= num_source
                    time=time+max(distance)/speed;%利用最大距离/速度得到截止到当前iteration的总时间
            end
            
           %% 在这里引入改进，禁区机制
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
                            %record_forbidden_center;
                   for k=1:size(record_forbidden_center,1)
                       forbiddenUb (k,:);
                       forbiddenLb(k,:);
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
                [fnew,best,nest,fitness]=get_best_nest2(nest,nest,fitness);%驱逐出去后比较一下找best
            end

            %% 找到截止目前的最佳位置bestnest和fmin函数
            if fnew<fmin
                fmin=fnew;
                bestnest=best;
            end
            tEnd=toc(tStart);
            
        end
            
                   if num_source == size(record_forbidden_center, 1)
                success = 1;
            else
                success = 0;
                   end

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
%             title( strcat('current generation=',num2str(T)), 'FontWeight', 'bold');
%             set(gca,'linewidth',0.5,'fontsize',8,'fontname','Times New Roman');
%              set(gcf,'position',[131 400 250 200])
% for i=1:size(record_forbidden_center,1)
%     hold on
%     center=record_forbidden_center(i,:);%中心
%     length=sqrt(2)/2*L; %边长
%     rotate=0*pi; %旋转角度
%     aa=[-1 -1 1 1 -1];
%     bb=[-1 1 1 -1 -1];
%     cc=complex(aa,bb)/sqrt(2);
%     cc=cc*exp(j*rotate)*length+complex(center(1),center(2));
%     line(real(cc),imag(cc),'LineWidth',1);
% end
% pause(0.1);
%     axis equal;
    end

