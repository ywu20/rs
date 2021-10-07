%% ��������������������ƵĸĽ���Ϊ�˷�ֹ����ֲ�����ֵ�㣬ͨ���ж�ÿһ��best��Χ�Ƿ��г���һ��������nest�����У����ⲿ��������Ϊ������
%% ��ֹ�Ժ�ɽ�ȥ��
%% Ԥ�ڽ���� �ҵ����оֲ���ֵ��


function [success time] = multi_source(initial_nest)
L=6;%���������α߳�

n_tolerance=5;%����n_tolerance��nest��best�������ж�Ϊ����
speed=3600;%ÿ���������˶��ٶ�m/min
time=0;%�����洢Max_condition�ֲ�ͬ�����N���ظ�ʵ���ʱ��
n=size(initial_nest, 1); %�������ĸ���
nd=2; %ÿ���񳲵�����ά�ȣ�������������ÿ���������Ķ�ά����x��y
pa=0.25; %�������ֲ����񵰵ĸ���
num_source = 4;
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
%% �����ռ�綨
M=50;%������������
Lb=1*ones(1,nd);  %x,y�������Сֵ
Ub=M*ones(1,nd);  %x,y��������ֵ
T=0;
        record_forbidden_center=[];%��¼������������\
        forbiddenUb=[];
forbiddenLb=[];
forbidden_center=[];
nest=initial_nest;
   %ÿ�����鴫������ʼ��һ�������Ʊ�������
        actual_nest=nest; %��¼����������ʵλ�ã�����֮�����ÿ���ƶ��ľ���distance
        %% �õ���ǰ����Ѵ�����
        fitness=10^10*ones(n,1);
        [fmin,bestnest,nest,fitness]=get_best_nest2(nest,nest,fitness); %�ȵõ���ʼʱ��bestnest��
        %% ��ʼѭ������
tStart=tic;
tEnd=0;
        while (tEnd<0.15)%�����ǲ������д��������Ͻ��ֹͣ
     T=T+1;
           %% ��һ�����£�levy�����´�ָ���õ�sensor��������������
            new_nest=get_cuckoos(nest,bestnest,Lb,Ub,record_forbidden_center,L); 
           %% �����һ�����µ�ʵ����ÿ��sensor���ƶ�����=levy ����ָ��λ��-��ǰλ��
            distance_first=sqrt(sum((new_nest-actual_nest).^2,2));%�����һ���ľ���
           %% ��������
            [fnew,best,nest,fitness]=get_best_nest2(nest,new_nest,fitness);
            actual_nest=new_nest;%����ʵ��λ��

            %% �ڶ������£�pa���ʷ��ֲ��ƶ�
            new_nest=empty_nests(nest,Lb,Ub,pa);%�����ڶ����ƶ�ָ��
            % �����ƶ�
            [fnew,best,nest,fitness]=get_best_nest2(nest,new_nest,fitness);

            %% ����ڶ��������Լ���һ������ʱ��
            distance_second=sqrt(sum((new_nest-actual_nest).^2,2));%����ڶ����ľ���
            actual_nest=new_nest;%����ʵ��λ��
            %sensor��ͬ���ģ����㵱ǰһ������룬����ô�ʱ��
            distance=distance_first+distance_second;%����Ӻ͵õ���ǰiteration�ľ���
            done = 0;
            if(size(record_forbidden_center,1)==num_source)
            for i=1:size(record_forbidden_center,1)
                if sum(abs(bestnest-[x0 y0]))<=5||sum(abs(bestnest-[x1 y1]))<= 5 ||sum(abs(bestnest-[x2 y2]))<= 5 || sum(abs(bestnest-[x3 y3]))<=5 || sum(abs(bestnest-[x4 y4]))<=5||sum(abs(bestnest-[x5 y5]))<= 5 ||sum(abs(bestnest-[x6 y6]))<= 5    
                   done=done+1;
                end
            end
            end
              if done ~= num_source
                    time=time+max(distance)/speed;%����������/�ٶȵõ���ֹ����ǰiteration����ʱ��
            end
            
           %% ����������Ľ�����������
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
                [fnew,best,nest,fitness]=get_best_nest2(nest,nest,fitness);%�����ȥ��Ƚ�һ����best
            end

            %% �ҵ���ֹĿǰ�����λ��bestnest��fmin����
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
%             title( strcat('current generation=',num2str(T)), 'FontWeight', 'bold');
%             set(gca,'linewidth',0.5,'fontsize',8,'fontname','Times New Roman');
%              set(gcf,'position',[131 400 250 200])
% for i=1:size(record_forbidden_center,1)
%     hold on
%     center=record_forbidden_center(i,:);%����
%     length=sqrt(2)/2*L; %�߳�
%     rotate=0*pi; %��ת�Ƕ�
%     aa=[-1 -1 1 1 -1];
%     bb=[-1 1 1 -1 -1];
%     cc=complex(aa,bb)/sqrt(2);
%     cc=cc*exp(j*rotate)*length+complex(center(1),center(2));
%     line(real(cc),imag(cc),'LineWidth',1);
% end
% pause(0.1);
%     axis equal;
    end

