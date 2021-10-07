%asynchronous cuckoo search for finding one absolute maximum
%function compare(initial_nest)
global c;
%% ÊµÑéµÄ²ÎÊıÉè¼Æ
speed=3600;%Ã¿¸ö´«¸ĞÆ÷ÔË¶¯ËÙ¶Èm/min
N=100;%Ã¿ÖÖÇé¿öÔËĞĞN´Î
time=zeros(N);%ÓÃÀ´´æ´¢Max_conditionÖÖ²»Í¬Çé¿öÏÂN´ÎÖØ¸´ÊµÑéµÄÊ±¼ä
successtime0 = 0; %¼ÇÂ¼³É¹¦´ÎÊı
successtime1=0;
L=6;%½ûÇøÕı·½ĞÎ±ß³¤
num_source = 4;
n_tolerance=5;%³¬¹ın_tolerance¸önestÔÚbest¸½½üÔòÅĞ¶ÏÎª½ûÇø
record_forbidden_center=[];
%% Ëã·¨µÄ²ÎÊıÉè¼Æ
n=10; %´«¸ĞÆ÷µÄ¸öÊı
nd=2; %Ã¿¸öÄñ³²µÄÃèÊöÎ¬¶È£¬¾ßÌåÎÊÌâÖĞÊÇÃ¿¸ö´«¸ĞÆ÷µÄ¶şÎ¬×ø±êx£¬y
pa=0.25; %ËŞÖ÷Äñ·¢ÏÖ²¼¹ÈÄñµ°µÄ¸ÅÂÊ
cuckooInfo = zeros(n,8);
Tmax=3000;
%% ÆøÌåÅ¨¶ÈµÄ²ÎÊıÉè¼Æ
k=0.03;v=0.03;
q0=(2*pi*k)^0.5;x0=10;y0=10; %Ä¿±êÔ´
q1=0.8*q0;      x1=23;y1=45; %¸ÉÈÅÔ´1
q2=0.7*q0;      x2=45;y2=10; %¸ÉÈÅÔ´2
q3=0.6*q0;      x3=40;y3=45; %¸ÉÈÅÔ´3
q4=0.5*q0;      x4=30;y4=30; %¸ÉÈÅÔ´4
q5=0.8*q0;      x5=13;y5=35; %¸ÉÈÅÔ´5
q6=0.7*q0;      x6=25;y6=16; %¸ÉÈÅÔ´6
q7=0.6*q0;      x7=14;y7=25; %¸ÉÈÅÔ´7
c=gas_concentration_generate(k,v,q0,x0,y0,q1,x1,y1,q2,x2,y2,q3,x3,y3,q4,x4,y4,q5,x5,y5,q6,x6,y6,q7,x7,y7);%¶àµãÔ´²úÉúµÄÅ¨¶È¿Õ¼ä

%% ËÑË÷¿Õ¼ä½ç¶¨
M=50;%ËÑË÷·½ÏòÉÏÏŞ
Lb=1*ones(1,nd);  %x,y·½ÏòµÄ×îĞ¡Öµ
Ub=M*ones(1,nd);  %x,y·½ÏòµÄ×î´óÖµ

%% far from source distribution
% for i=1:n
%     initial_nest(i,:)=4/5*Ub+Lb+(1/5*Ub-Lb).*rand(size(Lb)); %×¢£ºÕâÀïÍ¨¹ı³õÊ¼Î»ÖÃÉè¼Æ´«¸ĞÆ÷³õÊ¼Î»ÖÃ²»ºÃ£¬ÀëÔ´½ÏÔ¶µÄÇé¿ö
% end

%%uniform random distribution
% 
for i=1:n
initial_nest(i,:)=Lb+(Ub-Lb).*rand(size(Lb)); %×¢£ºÕâÀïÒ²¿ÉÒÔÍ¨¹ı³õÊ¼Î»ÖÃÉè¼Æ´«¸ĞÆ÷³õÊ¼Î»ÖÃ²»ºÃ£¬ÀëÔ´ºÜÔ¶µÄÇé¿ö 
end
successstat = zeros(1,2);
Max_condition=2;
time1 = zeros(1,N);
time2 = zeros(1,N);

 for m=1:Max_condition %²»Í¬Çé¿ö£¨ÆÕÍ¨csºÍ¸Ä½øcs£©
for j=1:N
   if m==1
       [successtime0 time1(1,j)] = main_2(initial_nest);
       if(successtime0 == 1)
           successstat(1,1) = successstat(1,1)+1;
        end 
       time1(1,j) = main_2(initial_nest);
   else
        [successtime1 time2(1,j)]=asynchronous_cuckoo_search(initial_nest);
   if(successtime1==1)
       successstat(1,2) = successstat(1,2) +1;
   end
        time2(1,j)=asynchronous_cuckoo_search(initial_nest);
   end
end
 end
% for j=1:N
%     if m==1
%       
%     time1(1,j)=main_2(initial_nest);
%        % successtime0=successtime0+b;
%         
%     else
%     time2(1,j)=asynchronous_cuckoo_search(initial_nest);
%     %successtime1=successtime1+b;
%     end
% end
% 
% end
% % % figure (2)
% % % boxplot (record_max_concentration)
% % for i=1:size(record_forbidden_center,1)
% %     hold on;
% %     center=record_forbidden_center(i,:);%ÖĞĞÄ
% %     length=sqrt(2)/2*L; %±ß³¤
% %     rotate=0*pi; %Ğı×ª½Ç¶È
% %     aa=[-1 -1 1 1 -1];
% %     bb=[-1 1 1 -1 -1];
% %     cc=complex(aa,bb)/sqrt(2);
% %     cc=cc*exp(j*rotate)*length+complex(center(1),center(2));
% %     line(real(cc),imag(cc),'LineWidth',1);
% % end
% % %     axis equal;
% %record_forbidden_center
% % Èı¸öÆÀ¼ÛÖ¸±ê£º
% % success_time %³É¹¦ÕÒµ½µÄ´ÎÊı¶Ô±È
% % sum(time(:,1)>time(:,2)) %N´ÎÊµÑéÖĞ£¬¼Ó°à¸Ä½øºóÓÃÊ±ÉÙµÄ´ÎÊı
% % average_time=mean(time) %N´ËÊµÑéÆ½¾ùÓÃÊ±¶Ô±
% % 
figure(4) 

  X = categorical({'classical','asynchronous'});
  successstat
 bar(X,successstat);
 title('Accuracy of Two Algorithms');
set(gca,'FontSize',14,'Fontname', 'Times New Roman');
% %  % %% »­Í¼
% % % 
% % 
figure (3)
TIME=zeros(N,Max_condition);
TIME(:,1)=time1;%ÇóN´ÎÖØ¸´ÊµÑéµÄÆ½¾ùÊ±¼ä
TIME(:,2) = time2;
% %plot(1:Max_condition,average_time,'.','MarkerSize',30)
boxplot (TIME, 'Labels', {'classical','asynchronous'});
set(gca,'FontSize',14,'Fontname', 'Times New Roman');
median(TIME(:,1))
median(TIME(:,2))
title('Searching Time');
%end