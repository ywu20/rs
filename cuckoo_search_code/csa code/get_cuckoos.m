function nest=get_cuckoos(nest,best,Lb,Ub,record_forbidden_center,L)

n=size(nest,1);
beta=3/2;
factor=10;
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
    nest(count,:)=simplebounds(s,Lb,Ub);
    count=count+1;
end
