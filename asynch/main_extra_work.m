%% ��������������������ƵĸĽ���Ϊ�˷�ֹ����ֲ�����ֵ�㣬ͨ���ж�ÿһ��best��Χ�Ƿ��г���һ��������nest�����У����ⲿ��������Ϊ������
%% ��ֹ�Ժ�ɽ�ȥ��
%% Ԥ�ڽ���� �ҵ����оֲ���ֵ��


clc;
clear;
global c;
%% ʵ��Ĳ������

Max_condition=3;%����Max_condition�ֲ�ͬ�������������׼�������֣���һ������ͨcs�㷨���ڶ�������������ĸĽ�
L=6;%���������α߳�
num_source = 4;
ratio=0.5;%�ж��ٱ���Ҫ�Ӱ�
n_tolerance=5;%����n_tolerance��nest��best�������ж�Ϊ����
speed=3600;%ÿ���������˶��ٶ�m/min
N=100;%ÿ���������N��
time=zeros(N,Max_condition);%�����洢Max_condition�ֲ�ͬ�����N���ظ�ʵ���ʱ��
successtime0 = 0;
successtime1 = 0;
successtime2=0;
% success_time=zeros(1,Max_condition);%������¼ÿһ�����������N�����м��γɹ��ҵ�
%% �㷨�Ĳ������
n=25; %�������ĸ���
nd=2; %ÿ���񳲵�����ά�ȣ�������������ÿ���������Ķ�ά����x��y
pa=0.25; %�������ֲ����񵰵ĸ���
%% ����Ũ�ȵĲ������
k=0.03;v=0.03;
q0=(2*pi*k)^0.5;x0=10;y0=10; %Ŀ��Դ
q1=0.8*q0;      x1=23;y1=45; %����Դ1
q2=0.7*q0;      x2=45;y2=10; %����Դ2
q3=0.6*q0;      x3=40;y3=45; %����Դ3
q4=0.5*q0;      x4=30;y4=30; %����Դ3
q5=0.8*q0;      x5=13;y5=35; %����Դ1
q6=0.7*q0;      x6=25;y6=16; %����Դ2
q7=0.6*q0;      x7=14;y7=25; %����Դ3
c=gas_concentration_generate(k,v,q0,x0,y0,q1,x1,y1,q2,x2,y2,q3,x3,y3,q4,x4,y4,q5,x5,y5,q6,x6,y6,q7,x7,y7);%���Դ������Ũ�ȿռ�
%% �㷨��ֹ����
Tmax=500; %������в���
record_max_concentration=zeros(Tmax,N);%������¼ÿ�����У�ÿ�ε��������Ũ�ȣ�����������ͼ��ok�����Ҫչ�ֲ�ͬ��pa��n����������԰ѱ�������Ϊ��ά��
%% �����ռ�綨
M=50;%������������
Lb=1*ones(1,nd);  %x,y�������Сֵ
Ub=M*ones(1,nd);  %x,y��������ֵ

%% ��������ʼλ�����
% for i=1:n
%     initial_nest(i,:)=4/5*Ub+Lb+(1/5*Ub-Lb).*rand(size(Lb)); %ע������ͨ����ʼλ����ƴ�������ʼλ�ò��ã���Դ��Զ�����
% end

%% ���������񳲣�������������ĳ�ʼλ�ã�Case 1��
for i=1:n
initial_nest(i,:)=Lb+(Ub-Lb).*rand(size(Lb)); %ע������Ҳ����ͨ����ʼλ����ƴ�������ʼλ�ò��ã���Դ��Զ�����
 end

%% ʵ�鿪ʼ
for m=1:Max_condition %��ͬ�������ͨcs�͸Ľ�cs��
    for j=1:N         %ÿ������ظ�����N��ʵ��
        record_forbidden_center=[];%��¼������������\
        forbiddenUb=[];
forbiddenLb=[];
forbidden_center=[];
nest=initial_nest;
   %ÿ�����鴫������ʼ��һ�������Ʊ�������
        actual_nest=nest; %��¼����������ʵλ�ã�����֮�����ÿ���ƶ��ľ���distance
        %% �õ���ǰ����Ѵ�����
        fitness=10^10*ones(n,1);
        [fmin,bestnest,nest,fitness]=get_best_nest(nest,nest,fitness); %�ȵõ���ʼʱ��bestnest��
        
        T=0;
        %% ��ʼѭ������
        
        while T<Tmax%�����ǲ������д��������Ͻ��ֹͣ
            T=T+1;
           %% ��һ�����£�levy�����´�ָ���õ�sensor��������������
             new_nest=get_cuckoos(nest,bestnest,Lb,Ub,record_forbidden_center,L);
           
           %% �����һ�����µ�ʵ����ÿ��sensor���ƶ�����=levy ����ָ��λ��-��ǰλ��
            distance_first=sqrt(sum((new_nest-actual_nest).^2,2));%�����һ���ľ���
           %% ��������
            [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness);
            actual_nest=new_nest;%����ʵ��λ��
            if m==3  %Ҳ���Լ���&T<Tmax/2
                new_nest=extra_work(new_nest,best,Lb,Ub,distance_first,speed,ratio);%ratio�������ȵ���ָ��λ�õ�sensor�ٽ���һ�����ƾ����levy����
            else
            %% �ڶ������£�pa���ʷ��ֲ��ƶ�
            new_nest=empty_nests(nest,Lb,Ub,pa);%�����ڶ����ƶ�ָ��
            end
            % �����ƶ�
            [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness);

            %% ����ڶ��������Լ���һ������ʱ��
            distance_second=sqrt(sum((new_nest-actual_nest).^2,2));%����ڶ����ľ���
            actual_nest=new_nest;%����ʵ��λ��
            %sensor��ͬ���ģ����㵱ǰһ������룬����ô�ʱ��
            distance=distance_first+distance_second;%����Ӻ͵õ���ǰiteration�ľ���
           if m==1 
               if max(record_max_concentration(:,j)) <= 4.5
            time(j,m)=time(j,m)+max(distance)/speed;%����������/�ٶȵõ���ֹ����ǰiteration����ʱ��
               else
                   time(j,m) = time(j,m);
               end
           end
           %% ����������Ľ�����������
           if m~=1
            if T>10 %ǰT/4�Ȳ����ǽ�������
                if fmin~=0
              cuckoo_index=cuckoo_division(nest,best,L);%����cuckoo_division��������ڵ�ǰbest��Χ��nest
                  if sum(cuckoo_index)>=n_tolerance %�������best��Χnest����n_tolerance������������
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
                           record_forbidden_center(end+1,:)=local_best;%��¼������������
                        forbiddenUb(end+1,:) = fUb;
                         forbiddenLb(end+1,:) = fLb;
                end%��¼������������
     end
 end
%                 break;
            
                    
                
              %% �ѽ����е�nest�߳�
                for num=1:n %ÿһ��nest����
                    for k=1:size(record_forbidden_center) %ÿһ����������
                     %% ����num��nest�ڵ�k�������ڣ��߳��������߳���Χ���Խ�������Ϊ���ģ��߳�Ϊ5L�������Σ�����������L�����Σ���Χ��
                        if forbiddenLb(k,1)<nest(num,1) & nest(num,1)<forbiddenUb(k,1) & forbiddenLb(k,2)<nest(num,2) & nest(num,2)<forbiddenUb(k,2)
%                             if nest(num,1)<record_forbidden_center(k,1) %��nest�ڽ�������
%                                 nest(num,1)=record_forbidden_center(k,1)-L/1-2*L*rand(1);%���nest������߳�
%                             else
%                                 nest(num,1)=record_forbidden_center(k,1)+L/2+2*L*rand(1);%���nest�����߳�
%                             end
nest(num,1)=record_forbidden_center(k,2)-5/2*L+5*L*rand(1);%���·���ҲҪ�仯������levy��������ѭ��
                            nest(num,2)=record_forbidden_center(k,2)-5/2*L+5*L*rand(1);%���·���ҲҪ�仯������levy��������ѭ��
                        end
                    end
                    nest(num,:)=simplebounds(nest(num,:),Lb,Ub);%Ҫ��֤�߳���nest��������Χ��
                    fitness(num)=fobj(nest(num,:));%����fitnessԭ��������fitness�ٺ�Ҳ��Ҫ��
                end
                [fnew,best,nest,fitness]=get_best_nest(nest,nest,fitness);%�����ȥ��Ƚ�һ����best
            end
if size (record_forbidden_center,1) ~= num_source 
    time(j,m)=time(j,m)+max(distance)/speed;
else
    time(j,m) = time(j,m);
end
           end
            %% �ҵ���ֹĿǰ�����λ��bestnest��fmin����
            if fnew<fmin,
                fmin=fnew;
                bestnest=best;
            end
            record_max_concentration(T,j)=-fmin;%�����洢ÿһ��ʵ��ÿһ��iteration�µ����Ũ��ֵ������֮��ͼ

%             %% ѭ����ֹ��������1��sensor��Ŀ��Դ����1m��Χ�ھ�ֹͣѭ��
%             if ((bestnest-[x0,y0]).^2)<1
%                 success_time(1,m)=success_time(1,m)+1;%�ɹ��ҵ����ɹ�����+1
%                 break;
%             end
%             figure(1);
%             hold off;
%             for i=1:n
%                 axis([0 M 0 M]);
%                 plot(nest(i,1),actual_nest(i,2),'*');%%n����ѡ�⣨nest���ò�ͬ��ɫ��*��ʾ
%                 plot(best(1,1),best(1,2),'.','markersize',20);%%���Ž⣨bestnest���üӴֵ�ʵ�ĵ��ʾ
%                 hold on;
%             end
%             x=1:M;
%             y=1:M;
%             global c;
%             contour(x,y,c);
%             for i=1:size(record_forbidden_center,1)
%                 hold on
%                 center=record_forbidden_center(i,:);%����
%                 length=sqrt(2)/2*L; %�߳�
%                 rotate=0*pi; %��ת�Ƕ�
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
          end% iteration����
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
     end     % �ظ�ʵ�����
end   % ��ͬ������� 
% 
% bestnest %��command�����չʾ����
% % figure (2)
% % boxplot (record_max_concentration)
% for i=1:size(record_forbidden_center,1)
%     hold on;
%     center=record_forbidden_center(i,:);%����
%     length=sqrt(2)/2*L; %�߳�
%     rotate=0*pi; %��ת�Ƕ�
%     aa=[-1 -1 1 1 -1];
%     bb=[-1 1 1 -1 -1];
%     cc=complex(aa,bb)/sqrt(2);
%     cc=cc*exp(j*rotate)*length+complex(center(1),center(2));
%     line(real(cc),imag(cc),'LineWidth',1);
% end
% %     axis equal;
record_forbidden_center
% ��������ָ�꣺
% success_time %�ɹ��ҵ��Ĵ����Ա�
% sum(time(:,1)>time(:,2)) %N��ʵ���У��Ӱ�Ľ�����ʱ�ٵĴ���
% average_time=mean(time) %N��ʵ��ƽ����ʱ�Ա� 
figure(4) 
successstat = [successtime0 successtime1 successtime2];
X=1:1:3;
bar(X,successstat);
title('Compare the number of times the requirement is met in 100 generations');
%% ��ͼ

figure (3)
%average_time=mean(time);%��N���ظ�ʵ���ƽ��ʱ��
%plot(1:Max_condition,average_time,'.','MarkerSize',30)
boxplot (time, 'Labels', {'conventional','improved', 'extra-improved'});
title('Compare the time it takes for the robots to meet the requirements');
