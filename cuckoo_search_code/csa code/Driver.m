
global c;
N=100;
time=zeros(N);
successtime0 = 0;
successtime1=0;
n=50;
nd=2;
M=50;
Lb=1*ones(1,nd); 
Ub=M*ones(1,nd);

k=0.03;v=0.03;
q0=(2*pi*k)^0.5;x0=10;y0=10;
q1=0.8*q0;      x1=23;y1=45;
q2=0.7*q0;      x2=45;y2=10;
q3=0.6*q0;      x3=40;y3=45;
q4=0.5*q0;      x4=30;y4=30;
q5=0.8*q0;      x5=13;y5=35;
q6=0.7*q0;      x6=25;y6=16;
q7=0.6*q0;      x7=14;y7=25;
c=gas_concentration_generate(k,v,q0,x0,y0,q1,x1,y1,q2,x2,y2,q3,x3,y3,q4,x4,y4,q5,x5,y5,q6,x6,y6,q7,x7,y7);

%% far from source distribution
% for i=1:n
%     initial_nest(i,:)=4/5*Ub+Lb+(1/5*Ub-Lb).*rand(size(Lb));
% end

%%uniform random distribution

for i=1:n
initial_nest(i,:)=Lb+(Ub-Lb).*rand(size(Lb));
end
successstat = zeros(1,2);
Max_condition=2;
time1 = zeros(1,N);
time2 = zeros(1,N);

 for m=1:Max_condition
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
  X = categorical({'classical','multi_source'});
  successstat
 bar(X,successstat);
 title('Accuracy of Two Algorithms');
set(gca,'FontSize',14,'Fontname', 'Times New Roman');

figure (3)
TIME=zeros(N,Max_condition);
TIME(:,1)=time1;
TIME(:,2) = time2;
boxplot (TIME, 'Labels', {'classical','multi-source'});
set(gca,'FontSize',14,'Fontname', 'Times New Roman');
median(TIME(:,1))
median(TIME(:,2))
title('Searching Time');

