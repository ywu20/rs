%% �������ܣ�Ϊ�˲��úܶ�nest����ֲ������У���ʵ�ֲ����񻮷�
%% �������ݣ���best�񳲱Ƚϣ�����bestΪ���ģ��߳�ΪL�����������г���n_tolerance��nest�Ļ����Ͱ��⼸��nest��ǣ������ָ�����Ϊ����
function cuckoo_index=cuckoo_division(nest,best,L)
N=size(nest,1);%nest����
cuckoo_index=zeros(N,1);%cuckoo_index��һ�������������ĵ�i�ж�Ӧ��i��nest�Ƿ��ڻ��ֵĽ����ڣ������㣬������һ
for i=1:N
    %% �����nest����bestΪ����LΪ�߳�����������
    if (best(1,1)-L/2)<nest(i,1) & nest(i,1)<(best(1,1)+L/2) & (best(1,2)-L/2)<nest(i,2) & nest(i,2)<(best(1,2)+L/2)
        cuckoo_index(i,1)=1;
    end
end