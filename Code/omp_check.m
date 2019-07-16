%checks the exact recovery condition
summ=0;
for i=1:length(dictionary(1,:))
    check=sum(abs(pinv(dictionary(:,initial_index_set))*dictionary(:,i)));
    if (check<1) && sum(i==index_set)==0
        summ=summ+1;
    end
end