clc
clear
%DCT dictionary creation(just the dictionary part found online)
bb=8;
eps=0.01;
Pn=ceil(sqrt(441));
DCT=zeros(bb,Pn);
for k=0:1:Pn-1
    V=cos([0:1:bb-1]'*k*pi/Pn);
    if k>0 
        V=V-mean(V); 
    end
    DCT(:,k+1)=V/norm(V);
end
DCT=kron(DCT,DCT);
%getting patches of the image
load('im_data_matrix.mat');
training_data_matrix=training_data_matrix(:,randperm(12474,11000));
sparsity_level=10;
representation_vectors=zeros(441,11000);
dictionary=DCT;
approximation=zeros(size(dictionary));
%first OMP
for i=1:11000
    i;
    %[approximation_data,estimate_ideal_signal,index_set,residual,chosen_atoms]=omp(dictionary,training_data_matrix(:,i),sparsity_level);
    [approximation(:,i),estimate_ideal_signal,index_set,residual]=omp_less_sparse(dictionary,training_data_matrix(:,i),sparsity_level);
    representation_vectors(:,i)=estimate_ideal_signal;
end
%
j=1;
say=0;
%fro_norm is frobenius norm
old_fro_norm=-1;
while(j<=180)
    for i=1:11000
        i;
        %[approximation_data,estimate_ideal_signal,index_set,residual,chosen_atoms]=omp(dictionary,training_data_matrix(:,i),sparsity_level);
        [approximation_data,estimate_ideal_signal,index_set,residual,chosen_atoms]=omp_less_sparse(dictionary,training_data_matrix(:,i),sparsity_level);
        %continue with old respresentation if OMP recovers a worse one
        if norm(training_data_matrix(:,i)-approximation_data)<norm(training_data_matrix(:,i)-approximation(:,i))
            say=say+1;
            representation_vectors(:,i)=estimate_ideal_signal;
        end
    end
    %not include DC
    for k=2:length(dictionary(1,:))
        error=training_data_matrix-dictionary*representation_vectors;
        indices_wk=find(representation_vectors(k,:)~=0);
        if length(indices_wk)>0
            shrinked_error_atom_moved=error(:,indices_wk)+dictionary(:,k)*representation_vectors(k,indices_wk);
            %minimization
            [U,S,V] = svds(shrinked_error_atom_moved,1,'largest');
            dictionary(:,k)=U;
            representation_vectors(k,indices_wk)=S*V';
            error(:,indices_wk)=shrinked_error_atom_moved-dictionary(:,k)*representation_vectors(k,indices_wk);
        %else
         %   unused_atom=unused_atom+1;
        end
    end
    approximation=dictionary*representation_vectors;
    fro_norm=norm((training_data_matrix-dictionary*representation_vectors),'fro');
    if(abs(fro_norm-old_fro_norm)<eps)
        break;
    end
    old_fro_norm=fro_norm;
    j=j+1;
end


