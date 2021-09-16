%classical CSA
function [success time] = main_2(initial_nest)
%there are some dummy variables just so that I can use same get_cuckoos for both algorithms.
L=6;
num_source = 4;
n_tolerance=3;
speed=3600;

n=size(initial_nest,1);

nd=2;
pa=0.25;
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

time =0;
M=50;
Lb=1*ones(1,nd);
Ub=M*ones(1,nd);
record_forbidden_center=[];

nest=initial_nest;
        actual_nest=nest;
        fitness=10^10*ones(n,1);
        [fmin,bestnest,nest,fitness]=get_best_nest2(nest,nest,fitness);
        T=0;
        tStart=tic;
        tEnd=0;
        while (tEnd<0.15) %tEnd is for without graphs, result may vary from computer from computer because every computer run with different speed. 
            new_nest=get_cuckoos(nest,bestnest,Lb,Ub,record_forbidden_center,L);
            distance_first=sqrt(sum((new_nest-actual_nest).^2,2));

            [fnew,bestnest,nest,fitness]=get_best_nest2(nest,new_nest,fitness);
            actual_nest=new_nest;

            new_nest=empty_nests(nest,Lb,Ub,pa);
            [fnew,bestnest,nest,fitness]=get_best_nest2(nest,new_nest,fitness);

            distance_second=sqrt(sum((new_nest-actual_nest).^2,2));
            actual_nest=new_nest;
            distance=distance_first+distance_second;
            if sum(abs(bestnest-[x0 y0]))>= 2
            time=time+max(distance)/speed;
            end
            tEnd=toc(tStart);
            
%             figure(1); %animates the procedure
%             hold off;
%             for i=1:n
%                 axis([0 M 0 M]);
%                 plot(nest(:,1),nest(:,2),'*');
%                 plot(bestnest(1,1),bestnest(1,2),'.','markersize',20);
%                 hold on;
%             end
%             x=1:M;
%             y=1:M;
%             contour(x,y,c);
%             pause(0.1);
        end

        if sum(abs(bestnest-[x0 y0]))<= 2
            success = 1;
        else
            success = 0;
        end
        end
