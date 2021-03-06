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
global c;
%% 搜索空间界定
M =50;
Lb=1*ones(1,nd);  %x,y方向的最小值
Ub=M*ones(1,nd);  %x,y方向的最大值

%% 传感器初始位置设计 Fig. 9,10
% for i=1:n
%     initial_nest(i,:)=4/5*Ub+Lb+(1/5*Ub-Lb).*rand(size(Lb)); %注：这里通过初始位置设计传感器初始位置不好，离源较远的情况
% end

%% 产生各个鸟巢（传感器）随机的初始位置（Case 1）Fig. 6,7
for i=1:size(initial_nest,1)
cuckooInfo(:,6) = fobj(initial_nest(i,:));
end

%% 实验开始
        cuckooInfo(:,2:3) = initial_nest;
        cuckooInfo(:,4:5) = initial_nest;
        %% 得到当前的最佳传感器
         [fmin,K]=min(cuckooInfo(:,6));
          bestnest= cuckooInfo(K,2:3);    %先得到初始时的bestnest等
        
%         %% 第一次飞行目标
%         %第一步更新：levy飞行下达指令（最好的sensor不动！！！！）
%         new_nest=get_cuckoos(cuckooInfo(:,2:3),bestnest,Lb,Ub,record_forbidden_center,L);
%         
%          % 第二步更新：pa概率发现并移动,记录目标位置
%          cuckooInfo(:,4:5)=empty_nests(new_nest,Lb,Ub,pa);%产生第二步移动指令
%          
%         %% 计算第一次飞行距离
%          distance_first=sqrt(sum((cuckooInfo(:,4:5)-cuckooInfo(:,2:3)).^2,2));
%          cuckooInfo(:,1)=distance_first;
% %         T=0;
% % while T<Tmax
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
%     while sum(abs(bestnest-[10 10])<= [2 2])~= 2%fobj(bestnest)>=-4 %没达到要求就继续找
  a= sqrt((bestnest(1,1)-10) * (bestnest(1,1)-10)) + ((bestnest(1,2)-10)* (bestnest(1,2)-10));
while a<= 1
[fmin first_arrive] =min(cuckooInfo(:,1)); %最先到达的鸟
            dest = get_cuckoos(cuckooInfo(:,2:3),bestnest,Lb,Ub,record_forbidden_center,L);
          
            dest = empty_nests(cuckooInfo(:,2:3), Lb, Ub, pa);
          
            cuckooInfo(first_arrive,4:5) = dest (first_arrive,:);
            d = sqrt(((cuckooInfo(first_arrive,2)-cuckooInfo(first_arrive,4))*(cuckooInfo(first_arrive,2)-cuckooInfo(first_arrive,4))+(cuckooInfo(first_arrive,3)-cuckooInfo(first_arrive,5))*(cuckooInfo(first_arrive,3)-cuckooInfo(first_arrive,5))));
            cuckooInfo(first_arrive,1) = cuckooInfo(first_arrive,1)+d;
            if fobj(cuckooInfo(first_arrive,4:5))<= fmin
                cuckooInfo(first_arrive,2:3) = cuckooInfo(first_arrive,4:5);
                cuckooInfo(first_arrive,6) = fobj(cuckooInfo(first_arrive,4:5));
                if cuckooInfo(first_arrive,6) <= fobj(bestnest)
                bestnest = cuckooInfo(first_arrive,2:3);
                end
            else
                cuckooInfo(first_arrive,1) = cuckooInfo(first_arrive,1)+d;
            end
            
%             cuckooInfo(first_arrive,2:3)=cuckooInfo(first_arrive,4:5); %飞到新位置
% 
%             %计算下一个目的地
%             if cuckooInfo(first_arrive,8) == 0
%                 [updated_first_arrive fitness] = get_best_nest(cuckooInfo(first_arrive,6:7),cuckooInfo(first_arrive,2:3),fitness);
%                 cuckooInfo (first_arrive,6:7)=updated_first_arrive; 
%                 cuckooInfo(first_arrive,8) = 1;
%             else
%                 updated_nest=get_cuckoos(cuckooInfo(:,6:7),bestnest,Lb,Ub,record_forbidden_center,L);
%                 updated_first_arrive = updated_nest(first_arrive,:);
%  %计算被不被扔掉              
%                 cuckoo = empty_nests(cuckooInfo(:,2:3),Lb,Ub,pa); %看新位置会不会比原来更好
%                 updated_first_arrive= cuckoo(first_arrive,:);
%                 cuckooInfo(first_arrive,8)= 0;
%             end
%             cuckooInfo(first_arrive,4:5) = updated_first_arrive;
% 
%            %如果新位置更好，变成best
%            if(fobj(cuckooInfo(first_arrive,2:3))<= fobj(bestnest))
%                bestnest=cuckooInfo(first_arrive,2:3);    
%            end
%            
%           %计算距离
%            distance_first_arrive = sqrt(sum((cuckooInfo(first_arrive,4:5)-cuckooInfo(first_arrive,2:3)).^2,2)); %计算新目标距离
%            cuckooInfo(first_arrive,1)= cuckooInfo(first_arrive,1) + distance_first_arrive;
           
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
%       
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
  a= sqrt((bestnest(1,1)-10) * (bestnest(1,1)-10)) + ((bestnest(1,2)-10)* (bestnest(1,2)-10));
    %if a(1,1)<= 2  && a(1,2)<= 2 
     %       success = 1;
      %      break;
       % else 
        %    success = 0;
     %end
tEnd=toc(tStart);

         end
   success = 1;
          
    time = max(cuckooInfo(:,1)')/speed;
   
       % tEnd
    % tEnd=toc(tStart);
    end