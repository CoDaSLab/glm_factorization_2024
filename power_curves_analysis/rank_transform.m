function Xtrans = rank_transform(X)

Xtrans = X;
for i=1:size(X,2)
    uX = unique(X(:,i));
    uX = sort(uX);
    jj=0;
    for j=1:length(uX)
        ind2 = find(X(:,i)==uX(j));
        Xtrans(ind2,i) = jj+mean(1:length(ind2));
        jj = jj+length(ind2);
    end
end