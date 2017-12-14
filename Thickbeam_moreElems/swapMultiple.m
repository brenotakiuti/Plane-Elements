function [y,order] = swapMultiple(y,order,ii,jj)

    swp = y([ii ii+1],:);
    y([ii ii+1],:) = y([jj jj+1],:) ;
    y([jj jj+1],:) = swp;
    swp = order([ii ii+1]);
    order([ii ii+1]) = order([jj jj+1]) ;
    order([jj jj+1]) = swp;
end