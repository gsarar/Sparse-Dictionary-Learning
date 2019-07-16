%load('globalTrainedDictionary.mat')
%dictionary=currDictionary;
%Creating DCT dictionary just the DCT dicitonary taken online
bb=8;
Pn=ceil(sqrt(441));
DCT=zeros(bb,Pn);
for k=0:1:Pn-1,
    V=cos([0:1:bb-1]'*k*pi/Pn);
    if k>0, V=V-mean(V); end;
    DCT(:,k+1)=V/norm(V);
end;
DCT=kron(DCT,DCT);
dictionary=DCT;

initial_index_set=randi(441,[1 3]);
data_vector=(sum(dictionary(:,initial_index_set)'))';
%beginning of algorithm
%initialization
sparsity_level=3;
[approximation_data,estimate_ideal_signal,index_set,residual]=omp(dictionary,data_vector,sparsity_level);
%
sum(abs(data_vector-approximation_data))
omp_check
summ