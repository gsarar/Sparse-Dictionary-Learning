function[approximation_data,estimate_ideal_signal,index_set,residual,chosen_atoms]=omp_error_bound(dictionary,data_vector,error)
residual=data_vector;
index_set=[];
t=1;
residual_norm_square=50;
while(residual_norm_square>error)
    max=0;
    index=0;
    for i=1:length(dictionary(1,:))
        val=abs(dot(residual,dictionary(:,i)));
        if val>max
            max=val;
            index=i;
        end
    end
    if index==0
        approximation_data=chosen_atoms*x_t;
        estimate_ideal_signal=zeros(length(dictionary(1,:)),1);
        estimate_ideal_signal(index_set)=x_t;
        break;
    end
    index_set=[index_set index];
    l=length(index_set);
    chosen_atoms=dictionary(:,index_set);
    x_t=pinv(chosen_atoms)*data_vector;
    residual=data_vector-chosen_atoms*x_t;
    residual_norm_square=(norm(residual)).^2;
    t=t+1;
end
approximation_data=chosen_atoms*x_t;
estimate_ideal_signal=zeros(length(dictionary(1,:)),1);
estimate_ideal_signal(index_set)=x_t;
end
