%asynchronous cuckoo search for finding one absolute maximum
function [success time] = asynchronous_cuckoo_search(initial_nest) 
%% 实验的参数设计
speed=3600;%每个传感器运动速度m/min
N=size(initial_nest,1);%每种情况运行N次
time=zeros(N);%用来存储Max_condition种不同情况下N次重复实验的时间
successtime0 = 0; %记录成功次数
L=6;%禁区正方形边长
num_source = 4;
n_tolerance=5;%超过n_tolerance个nest在best附近则判断为禁区
record_forbidden_center=[];
n=size(initial_nest,1);
%% 算法的参数设计
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

%% 传感器初始位置设计 Fig. 9,10
% for i=1:n
%     initial_nest(i,:)=4/5*Ub+Lb+(1/5*Ub-Lb).*rand(size(Lb)); %注：这里通过初始位置设计传感器初始位置不好，离源较远的情况
% end

%% 产生各个鸟巢（传感器）随机的初始位置（Case 1）Fig. 6,7
for i=1:size(initial_nest,1)
fitness(i) = fobj(initial_nest(i,:));
end

%% 实验开始
        cuckooInfo(:,2:3) = initial_nest;
        cuckooInfo(:,6:7) = initial_nest;
        %% 得到当前的最佳传感器
         [fmin,K]=min(fitness);
          bestnest=cuckooInfo(K,2:3);    %先得到初始时的bestnest等
        
        %% 第一次飞行目标
        %第一步更新：levy飞行下达指令（最好的sensor不动！！！！）
        new_nest=get_cuckoos(cuckooInfo(:,2:3),bestnest,Lb,Ub,record_forbidden_center,L);
        
         % 第二步更新：pa概率发现并移动,记录目标位置
         cuckooInfo(:,4:5)=empty_nests(new_nest,Lb,Ub,pa);%产生第二步移动指令
         
        %% 计算第一次飞行距离
         distance_first=sqrt(sum((cuckooInfo(:,4:5)-cuckooInfo(:,2:3)).^2,2));
         cuckooInfo(:,1)=distance_first;
%         T=0;
% while T<Tmax
%     T=T+1;
%abs(bestnest-[10 10])<= [2 2]
  tStart=tic;
        tEnd=0;
        bestnest = cuckooInfo(1,2:3);
%
%while tEnd<0.15 %* 30 %(tEnd<0.1* size(initial_nest,1))
        T=0;
     %while T < 500 * size(initial_nest)
  
            T=T+1;
     while sum(abs(bestnest-[10 10])<= [2 2])~= 2%fobj(bestnest)>=-4 %没达到要求就继续找
            [fmin first_arrive] =min(cuckooInfo(:,1)); %最先到达的鸟
            cuckooInfo(first_arrive,2:3)=cuckooInfo(first_arrive,4:5); %飞到新位置

            %计算下一个目的地
            if cuckooInfo(first_arrive,8) == 0
                [updated_first_arrive fitness] = get_best_nest(cuckooInfo(first_arrive,6:7),cuckooInfo(first_arrive,2:3),fitness);
                cuckooInfo (first_arrive,6:7)=updated_first_arrive; 
                cuckooInfo(first_arrive,8) = 1;
            else
                updated_nest=get_cuckoos(cuckooInfo(:,6:7),bestnest,Lb,Ub,record_forbidden_center,L);
                updated_first_arrive = updated_nest(first_arrive,:);
 %计算被不被扔掉              
                cuckoo = empty_nests(cuckooInfo(:,2:3),Lb,Ub,pa); %看新位置会不会比原来更好
                updated_first_arrive= cuckoo(first_arrive,:);
                cuckooInfo(first_arrive,8)= 0;
            end
            cuckooInfo(first_arrive,4:5) = updated_first_arrive;

           %如果新位置更好，变成best
           if(fobj(cuckooInfo(first_arrive,2:3))<= fobj(bestnest))
               bestnest=cuckooInfo(first_arrive,2:3);    
           end
           
          %计算距离
           distance_first_arrive = sqrt(sum((cuckooInfo(first_arrive,4:5)-cuckooInfo(first_arrive,2:3)).^2,2)); %计算新目标距离
           cuckooInfo(first_arrive,1)= cuckooInfo(first_arrive,1) + distance_first_arrive;
           
    %       if fobj(bestnest)>= 4.5 
  %successtime0 = successtime0 + 1;
 % break;
    %       end
%             figure(1);
%             hold off;
%             for i=1:n
%                 axis([0 M 0 M]);
%                 plot(cuckooInfo(:,2),cuckooInfo(:,3),'*');%%n个备选解（nest）用不同颜色的*表示
%                 plot(bestnest(1,1),bestnest(1,2),'.','markersize',20);%%最优解（bestnest）用加粗的实心点表示
%                 hold on;
%             end
%             x=1:M;
%             y=1:M;
%             global c;
%             contour(x,y,c);
%             for i=1:size(record_forbidden_center,1)
%                 hold on
%                 center=record_forbidden_center(i,:);%中心
%                 length=sqrt(2)/2*L; %边长2e546t4
%                 rotate=0*pi; %旋转角度
%                 aa=[-1 -1 1 1 -1];
%                 bb=[-1 1 1 -1 -1];
%                 cc=complex(aa,bb)/sqrt(2);
%                 cc=cc*exp(j*rotate)*length+complex(center(1),center(2));
%                 line(real(cc),imag(cc),'LineWidth',1);
%             end
%            title(num2str(T));
%           pause(0.1);
tEnd=toc(tStart);

         end

          
    time = max(cuckooInfo(:,1)')/speed;
    a= abs(bestnest-[10 10]);
    if a(1,1)<= 2  && a(1,2)<= 2 
            success = 1;
        else 
            success = 0;
     end
       % tEnd
    % tEnd=toc(tStart);
    end