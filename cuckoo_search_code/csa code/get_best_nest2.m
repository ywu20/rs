function [fmin best nest fitness]=get_best_nest(nest,newnest,fitness)

for j=1:size(nest,1),
    fnew=fobj(newnest(j,:));
    if fnew<=fitness(j) 
    fitness(j) = fnew;
    nest(j,:)=newnest(j,:);
    end
end

[fmin,K]=min(fitness) ;
best=nest(K,:);