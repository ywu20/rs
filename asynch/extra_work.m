function nest=extra_work(nest,best,Lb,Ub,distance,speed,ratio)

n=floor(size(nest,1)*ratio);%n����Ҫ�Ӱ��nest������ȡ��
%% ȷ����Ҫ�Ӱ��nest��nest�����е�λ��
sorted_distance=sort(distance);%�����д�С��������
% m=find(distance<sorted_distance(n))
k=zeros(size(nest,1),1);%Ϊ�˵õ�һ����0,1��ɵ������������üӰ��nest��Ӧ0����Ҫ�Ӱ�Ķ�Ӧ1
for i=1:size(nest,1)
    if distance(i)<=sorted_distance(n)
        k(i)=1;
    end
end
%% ����metagna�㷨���levy���мӰ�
beta=3/2;
factor=1;%%%������Ե�����С
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
%% ȷ��ʵ�ʼӰ��nestλ��
for i=1:size(nest,1)
    s=nest(i,:);
    u=randn(size(s))*sigma;
    v=randn(size(s));
    step=u./abs(v).^(1/beta);
    stepsize=factor*step.*(s-best);
    s=s+stepsize.*randn(size(s))*k(i);%���Ӱ��k��i��=0������levy���мӰ�
    nest(i,:)=simplebounds(s,Lb,Ub);
end        
