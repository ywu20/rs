function nest=get_cuckoos(nest,best,Lb,Ub,record_forbidden_center,L)
% Levy flights
n=size(nest,1);
% Levy exponent and coefficient
% For details, see equation (2.21), Page 16 (chapter 2) of the book
% X. S. Yang, Nature-Inspired Metaheuristic Algorithms, 2nd Edition, Luniver Press, (2010).
beta=3/2;
factor=10;%%%这里可以调整大小
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
count=1;
while count<=n
    flag=0;
    s=nest(count,:);
    u=randn(size(s))*sigma;
    v=randn(size(s));
    step=u./abs(v).^(1/beta);
    stepsize=factor*step.*(s-0.9999999.*best);
    s=s+stepsize.*randn(size(s));
%     for i=1:size(record_forbidden_center,1)
%         if (record_forbidden_center(i,1)-L/2)<s(1,1) & s(1,1)<(record_forbidden_center(i,1)+L/2) & (record_forbidden_center(i,2)-L/2)<s(1,2) & s(1,2)<(record_forbidden_center(i,2)+L/2)
%             flag=1;
%             break;
%         end
%     end
%     if flag==1
%         continue%若levy产生的指令位置在任何禁区内，不执行下面两句回到while开始出再次产生新的指令
%     end
    nest(count,:)=simplebounds(s,Lb,Ub);
    count=count+1;
end