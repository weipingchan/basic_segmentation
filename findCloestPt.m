function cloestPt=findCloestPt(ptList,tarPt)
    %compute Euclidean distances:
    distances = sqrt(sum(bsxfun(@minus, ptList,tarPt).^2,2));
    %find the smallest distance and use that as an index into B:
    cloestPt0=ptList(distances==min(distances),:);

    if size(cloestPt0,1)>1
        cloestPt=cloestPt0(1,:);
    else
        cloestPt=cloestPt0;
    end
end