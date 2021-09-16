
function cuckoo_index=cuckoo_division(nest,best,L)
N=size(nest,1);
cuckoo_index=zeros(N,1);
for i=1:N
    %%
    if (best(1,1)-L/2)<nest(i,1) & nest(i,1)<(best(1,1)+L/2) & (best(1,2)-L/2)<nest(i,2) & nest(i,2)<(best(1,2)+L/2)
        cuckoo_index(i,1)=1;
    end
end
