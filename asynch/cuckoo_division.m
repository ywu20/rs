%% 函数功能：为了不让很多nest陷入局部最优中，先实现布谷鸟划分
%% 划分依据：和best鸟巢比较，若以best为中心，边长为L的正方形中有超过n_tolerance个nest的话，就把这几个nest标记，并划分该区域为禁区
function cuckoo_index=cuckoo_division(nest,best,L)
N=size(nest,1);%nest个数
cuckoo_index=zeros(N,1);%cuckoo_index是一个列向量，它的第i行对应第i个nest是否在划分的禁区内，是置零，不是置一
for i=1:N
    %% 如果有nest在以best为中心L为边长的正方形内
    if (best(1,1)-L/2)<nest(i,1) & nest(i,1)<(best(1,1)+L/2) & (best(1,2)-L/2)<nest(i,2) & nest(i,2)<(best(1,2)+L/2)
        cuckoo_index(i,1)=1;
    end
end