%asynchronous cuckoo search for finding one absolute maximum
function [success time] = asynchronous_cuckoo_search(initial_nest) 
%% ʵ��Ĳ������
speed=3600;%ÿ���������˶��ٶ�m/min
N=size(initial_nest,1);%ÿ���������N��
time=zeros(N);%�����洢Max_condition�ֲ�ͬ�����N���ظ�ʵ���ʱ��
successtime0 = 0; %��¼�ɹ�����
L=6;%���������α߳�
num_source = 4;
n_tolerance=5;%����n_tolerance��nest��best�������ж�Ϊ����
record_forbidden_center=[];
n=size(initial_nest,1);
%% �㷨�Ĳ������
nd=2; %ÿ���񳲵�����ά�ȣ�������������ÿ���������Ķ�ά����x��y
pa=0.25; %�������ֲ����񵰵ĸ���
cuckooInfo = zeros(n,8);
Tmax=3000;
%% ����Ũ�ȵĲ������
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

%% �����ռ�綨
M=50;%������������
Lb=1*ones(1,nd);  %x,y�������Сֵ
Ub=M*ones(1,nd);  %x,y��������ֵ

%% ��������ʼλ����� Fig. 9,10
% for i=1:n
%     initial_nest(i,:)=4/5*Ub+Lb+(1/5*Ub-Lb).*rand(size(Lb)); %ע������ͨ����ʼλ����ƴ�������ʼλ�ò��ã���Դ��Զ�����
% end

%% ���������񳲣�������������ĳ�ʼλ�ã�Case 1��Fig. 6,7
for i=1:size(initial_nest,1)
fitness(i) = fobj(initial_nest(i,:));
end

%% ʵ�鿪ʼ
        cuckooInfo(:,2:3) = initial_nest;
        cuckooInfo(:,6:7) = initial_nest;
        %% �õ���ǰ����Ѵ�����
         [fmin,K]=min(fitness);
          bestnest=cuckooInfo(K,2:3);    %�ȵõ���ʼʱ��bestnest��
        
        %% ��һ�η���Ŀ��
        %��һ�����£�levy�����´�ָ���õ�sensor��������������
        new_nest=get_cuckoos(cuckooInfo(:,2:3),bestnest,Lb,Ub,record_forbidden_center,L);
        
         % �ڶ������£�pa���ʷ��ֲ��ƶ�,��¼Ŀ��λ��
         cuckooInfo(:,4:5)=empty_nests(new_nest,Lb,Ub,pa);%�����ڶ����ƶ�ָ��
         
        %% �����һ�η��о���
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
     while sum(abs(bestnest-[10 10])<= [2 2])~= 2%fobj(bestnest)>=-4 %û�ﵽҪ��ͼ�����
            [fmin first_arrive] =min(cuckooInfo(:,1)); %���ȵ������
            cuckooInfo(first_arrive,2:3)=cuckooInfo(first_arrive,4:5); %�ɵ���λ��

            %������һ��Ŀ�ĵ�
            if cuckooInfo(first_arrive,8) == 0
                [updated_first_arrive fitness] = get_best_nest(cuckooInfo(first_arrive,6:7),cuckooInfo(first_arrive,2:3),fitness);
                cuckooInfo (first_arrive,6:7)=updated_first_arrive; 
                cuckooInfo(first_arrive,8) = 1;
            else
                updated_nest=get_cuckoos(cuckooInfo(:,6:7),bestnest,Lb,Ub,record_forbidden_center,L);
                updated_first_arrive = updated_nest(first_arrive,:);
 %���㱻�����ӵ�              
                cuckoo = empty_nests(cuckooInfo(:,2:3),Lb,Ub,pa); %����λ�û᲻���ԭ������
                updated_first_arrive= cuckoo(first_arrive,:);
                cuckooInfo(first_arrive,8)= 0;
            end
            cuckooInfo(first_arrive,4:5) = updated_first_arrive;

           %�����λ�ø��ã����best
           if(fobj(cuckooInfo(first_arrive,2:3))<= fobj(bestnest))
               bestnest=cuckooInfo(first_arrive,2:3);    
           end
           
          %�������
           distance_first_arrive = sqrt(sum((cuckooInfo(first_arrive,4:5)-cuckooInfo(first_arrive,2:3)).^2,2)); %������Ŀ�����
           cuckooInfo(first_arrive,1)= cuckooInfo(first_arrive,1) + distance_first_arrive;
           
    %       if fobj(bestnest)>= 4.5 
  %successtime0 = successtime0 + 1;
 % break;
    %       end
%             figure(1);
%             hold off;
%             for i=1:n
%                 axis([0 M 0 M]);
%                 plot(cuckooInfo(:,2),cuckooInfo(:,3),'*');%%n����ѡ�⣨nest���ò�ͬ��ɫ��*��ʾ
%                 plot(bestnest(1,1),bestnest(1,2),'.','markersize',20);%%���Ž⣨bestnest���üӴֵ�ʵ�ĵ��ʾ
%                 hold on;
%             end
%             x=1:M;
%             y=1:M;
%             global c;
%             contour(x,y,c);
%             for i=1:size(record_forbidden_center,1)
%                 hold on
%                 center=record_forbidden_center(i,:);%����
%                 length=sqrt(2)/2*L; %�߳�2e546t4
%                 rotate=0*pi; %��ת�Ƕ�
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