function[approximation_data,estimate_ideal_signal,index_set,residual,chosen_atoms]=omp(dictionary,data_vector,sparsity_level)
residual=data_vector;
index_set=[];
t=1;
while(t<=sparsity_level)
    max=0;
    index=0;
    for i=1:length(dictionary(1,:))
        %finding closest atom
        val=abs(dot(residual,dictionary(:,i)));
        if val>max
            max=val;
            index=i;
        end
    end
    %gathering chosen atom indexes and finding chosen atoms
    index_set=[index_set index];
    chosen_atoms=dictionary(:,index_set);
    x_t=pinv(chosen_atoms)*data_vector;
    %finding residual
    residual=data_vector-chosen_atoms*x_t;
    t=t+1;
end
approximation_data=chosen_atoms*x_t;
estimate_ideal_signal=zeros(length(dictionary(1,:)),1);
%finding sparse representation
estimate_ideal_signal(index_set)=x_t;
end
