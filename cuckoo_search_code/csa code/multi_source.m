function [success time] = multi_source(initial_nest)
L=6;

n_tolerance=5;
speed=3600;
time=0;
n=size(initial_nest, 1); 
nd=2; 
pa=0.25;
num_source = 4; % change num_source settings here
k=0.03;v=0.03;
q0=(2*pi*k)^0.5;x0=10;y0=10; 
q1=0.8*q0;      x1=23;y1=45;
q2=0.7*q0;      x2=45;y2=10; 
q3=0.6*q0;      x3=40;y3=45; 
q4=0.5*q0;      x4=30;y4=30; 
q5=0.8*q0;      x5=13;y5=35; 
q6=0.7*q0;      x6=25;y6=16; 
q7=0.6*q0;      x7=14;y7=25; 
c=gas_concentration_generate(k,v,q0,x0,y0,q1,x1,y1,q2,x2,y2,q3,x3,y3,q4,x4,y4,q5,x5,y5,q6,x6,y6,q7,x7,y7);%多点源产生的浓度空间

M=50;
Lb=1*ones(1,nd);  
Ub=M*ones(1,nd);  
T=0;
record_forbidden_center=[];
forbiddenUb=[];
forbiddenLb=[];
forbidden_center=[];
nest=initial_nest;

actual_nest=nest; 

fitness=10^10*ones(n,1);
[fmin,bestnest,nest,fitness]=get_best_nest2(nest,nest,fitness); 
tStart=tic;
tEnd=0;
    while (tEnd<0.15) % this time is for without graph 
     T=T+1;
            new_nest=get_cuckoos(nest,bestnest,Lb,Ub,record_forbidden_center,L); 
            distance_first=sqrt(sum((new_nest-actual_nest).^2,2));
            [fnew,best,nest,fitness]=get_best_nest2(nest,new_nest,fitness);
            actual_nest=new_nest;

            new_nest=empty_nests(nest,Lb,Ub,pa);
            [fnew,best,nest,fitness]=get_best_nest2(nest,new_nest,fitness);
            distance_second=sqrt(sum((new_nest-actual_nest).^2,2));
            actual_nest=new_nest;
            distance=distance_first+distance_second;
            done = 0;
            if(size(record_forbidden_center,1)==num_source)
            for i=1:size(record_forbidden_center,1)
                if sum(abs(bestnest-[x0 y0]))<=5||sum(abs(bestnest-[x1 y1]))<= 5 ||sum(abs(bestnest-[x2 y2]))<= 5 || sum(abs(bestnest-[x3 y3]))<=5 || sum(abs(bestnest-[x4 y4]))<=5||sum(abs(bestnest-[x5 y5]))<= 5 ||sum(abs(bestnest-[x6 y6]))<= 5    
                   done=done+1;
                end
            end
            end
              if done ~= num_source
                    time=time+max(distance)/speed;
            end
            
            if T>10 
                if fmin~=0
              cuckoo_index=cuckoo_division(nest,best,L);
                  if sum(cuckoo_index)>=n_tolerance
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
                           record_forbidden_center(end+1,:)=local_best;
                        forbiddenUb(end+1,:) = fUb;
                         forbiddenLb(end+1,:) = fLb;
                   end
     end
                end
                for num=1:n 
                    for k=1:size(record_forbidden_center) 
                        if forbiddenLb(k,1)<nest(num,1) & nest(num,1)<forbiddenUb(k,1) & forbiddenLb(k,2)<nest(num,2) & nest(num,2)<forbiddenUb(k,2)
                            nest(num,1)=record_forbidden_center(k,2)-5/2*L+5*L*rand(1);
                            nest(num,2)=record_forbidden_center(k,2)-5/2*L+5*L*rand(1);
                        end
                    end
                    nest(num,:)=simplebounds(nest(num,:),Lb,Ub);
                    fitness(num)=fobj(nest(num,:));
                end
                [fnew,best,nest,fitness]=get_best_nest2(nest,nest,fitness);
            end

            if fnew<fmin
                fmin=fnew;
                bestnest=best;
            end
            tEnd=toc(tStart);

%            figure(2);
%             hold off;
%             for i=1:n
%                 axis([0 M 0 M]);
%                 plot(nest(i,1),actual_nest(i,2),'*');
%                 plot(best(1,1),best(1,2),'.','markersize',20);
%                 hold on;
%             end
%             x=1:M;
%             y=1:M;
%             contour(x,y,c);
%             for i=1:size(record_forbidden_center,1)
%                 hold on
%                 center=record_forbidden_center(i,:);
%                 length=sqrt(2)/2*L; 
%                 rotate=0*pi; 
%                 aa=[-1 -1 1 1 -1];
%                 bb=[-1 1 1 -1 -1];
%                 cc=complex(aa,bb)/sqrt(2);
%                 cc=cc*exp(j*rotate)*length+complex(center(1),center(2));
%                 line(real(cc),imag(cc),'LineWidth',1);
%             end
% pause(0.1);
    end
if num_source == size(record_forbidden_center, 1)
                success = 1;
            else
                success = 0;
                   end
    end