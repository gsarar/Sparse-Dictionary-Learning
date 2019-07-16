clc
clear
success_matrix=zeros(4,50);
time=0;
for n=[10 20 30 300]
    %for n=20
    time=time+1;
    for trial=1:50
        %for trial=1
        unused_atom=0;
        %generating dictionary
        generating_dictionary=rand([20 50]);
        %generating_dictionary=2.*rand([20 50])-1;
        %each column normalized l2 norm
        for i=1:length(generating_dictionary(1,:))
            generating_dictionary(:,i)=generating_dictionary(:,i)/norm(generating_dictionary(:,i));
        end
        %data signals
        data_signal=zeros(20,1500);
        data_signal_noise=zeros(20,1500);
        for k=1:1500
            coefficients=randperm(50,3);
            data_signal(:,k)=generating_dictionary(:,coefficients)*rand([3 1]);
            data_signal_noise(:,k)=awgn(data_signal(:,k),n,'measured');
        end
        %initializing dictionary
        initial_dictionary=data_signal(:,randperm(1500,50));
        for i=1:length(initial_dictionary(1,:))
            initial_dictionary(:,i)=initial_dictionary(:,i)/norm(initial_dictionary(:,i));
        end
        sparsity_level=3;
        representation_vectors=zeros(50,1500);
        dictionary=initial_dictionary;
        approximation=zeros(size(data_signal));
        %first OMP
        for i=1:1500
            i;
            [approximation(:,i),estimate_ideal_signal,index_set,residual]=omp_less_sparse(dictionary,data_signal_noise(:,i),sparsity_level);
            representation_vectors(:,i)=estimate_ideal_signal;
        end
        %
        j=1;
        say=0;
        while(j<=80)
            %OMP
            for i=1:1500
                i;
                [approximation_data,estimate_ideal_signal,index_set,residual]=omp_less_sparse(dictionary,data_signal_noise(:,i),sparsity_level);
                  %continue with old respresentation if OMP recovers a worse one
                if norm(data_signal_noise(:,i)-approximation_data)<norm(data_signal_noise(:,i)-approximation(:,i))
                    say=say+1;
                    representation_vectors(:,i)=estimate_ideal_signal;
                end
            end
            
            for k=1:length(dictionary(1,:))
                %%K-svd algortihm
                error=data_signal_noise-dictionary*representation_vectors;
                indices_wk=find(representation_vectors(k,:)~=0);
                if length(indices_wk)>0
                    shrinked_error_atom_moved=error(:,indices_wk)+dictionary(:,k)*representation_vectors(k,indices_wk);
                    %minimization-getting the U,V vectors corresponding to
                    %largest singular value and updating
                    [U,S,V] = svds(shrinked_error_atom_moved,1,'largest');
                    dictionary(:,k)=U;
                    representation_vectors(k,indices_wk)=S*V';
                    error(:,indices_wk)=shrinked_error_atom_moved-dictionary(:,k)*representation_vectors(k,indices_wk);
                else
                    unused_atom=unused_atom+1;
                end
            end
            approximation=dictionary*representation_vectors;
            j=j+1;
        end
        %evaluation
        success=0;
        for i=1:50
            [val,ind]=min(sum(abs(repmat(generating_dictionary(:,i),[1,50])-dictionary)));
            1-abs(generating_dictionary(:,i)'*dictionary(:,ind))
            %if (0.01>1-abs(generating_dictionary(:,i)'*dictionary(:,ind)))
            if (0.1>1-abs(generating_dictionary(:,i)'*dictionary(:,ind)))
                success=success+1;
            end
        end
        success_matrix(time,trial)=success;
    end
end