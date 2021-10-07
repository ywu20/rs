function nest=extra_work(nest,best,Lb,Ub,distance,speed,ratio)

n=floor(size(nest,1)*ratio);%n是需要加班的nest个数，取整
%% 确定需要加班的nest在nest变量中的位置
sorted_distance=sort(distance);%按照列从小到大排序
% m=find(distance<sorted_distance(n))
k=zeros(size(nest,1),1);%为了得到一个由0,1组成的行向量，不用加班的nest对应0，需要加班的对应1
for i=1:size(nest,1)
    if distance(i)<=sorted_distance(n)
        k(i)=1;
    end
end
%% 利用metagna算法设计levy飞行加班
beta=3/2;
factor=1;%%%这里可以调整大小
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
%% 确定实际加班后nest位置
for i=1:size(nest,1)
    s=nest(i,:);
    u=randn(size(s))*sigma;
    v=randn(size(s));
    step=u./abs(v).^(1/beta);
    stepsize=factor*step.*(s-best);
    s=s+stepsize.*randn(size(s))*k(i);%不加班的k（i）=0，不用levy飞行加班
    nest(i,:)=simplebounds(s,Lb,Ub);
end        
