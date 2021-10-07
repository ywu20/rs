%asynchronous cuckoo search for finding one absolute maximum
%function compare(initial_nest)
global c;
%% ʵ��Ĳ������
speed=3600;%ÿ���������˶��ٶ�m/min
N=100;%ÿ���������N��
time=zeros(N);%�����洢Max_condition�ֲ�ͬ�����N���ظ�ʵ���ʱ��
successtime0 = 0; %��¼�ɹ�����
successtime1=0;
L=6;%���������α߳�
num_source = 4;
n_tolerance=5;%����n_tolerance��nest��best�������ж�Ϊ����
record_forbidden_center=[];
%% �㷨�Ĳ������
n=50; %�������ĸ���
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

%% far from source distribution
% for i=1:n
%     initial_nest(i,:)=4/5*Ub+Lb+(1/5*Ub-Lb).*rand(size(Lb)); %ע������ͨ����ʼλ����ƴ�������ʼλ�ò��ã���Դ��Զ�����
% end

%%uniform random distribution
% 
for i=1:n
initial_nest(i,:)=Lb+(Ub-Lb).*rand(size(Lb)); %ע������Ҳ����ͨ����ʼλ����ƴ�������ʼλ�ò��ã���Դ��Զ����� 
end
successstat = zeros(1,2);
Max_condition=2;
time1 = zeros(1,N);
time2 = zeros(1,N);

 for m=1:Max_condition %��ͬ�������ͨcs�͸Ľ�cs��
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
% %  % %% ��ͼ
% % % 
% % 
figure (3)
TIME=zeros(N,Max_condition);
TIME(:,1)=time1;%��N���ظ�ʵ���ƽ��ʱ��
TIME(:,2) = time2;
% %plot(1:Max_condition,average_time,'.','MarkerSize',30)
boxplot (TIME, 'Labels', {'classical','multiSource'});
set(gca,'FontSize',14,'Fontname', 'Times New Roman');
median(TIME(:,1))
median(TIME(:,2))
title('Searching Time');
%end