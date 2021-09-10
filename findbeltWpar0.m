function areaBoundaryLength=findbeltWpar0(ferqSeries0,RF,forehindCornerA,realCen)
if RF=='R'
%For Right
    Rseries=ferqSeries0(ferqSeries0(:,1)>0,:);
    Rseries(:,1)=abs(Rseries(:,1));
    Rseries=sortrows(Rseries);
    panelLimit=floor(abs(forehindCornerA(1)-round(realCen(1))));
else
%For Left
    Rseries=ferqSeries0(ferqSeries0(:,1)<0,:);
    Rseries(:,1)=abs(Rseries(:,1));
    Rseries=sortrows(Rseries);
    panelLimit=floor(abs(forehindCornerA(1)-round(realCen(1))));
end

areaSumDiff=diff(Rseries(1:panelLimit,2));
areaSumDiffSmo=smoothdata(areaSumDiff,'loess');
areaSumDiff0=areaSumDiffSmo(1:end-1);
areaSumDiff1=areaSumDiffSmo(2:end);
areaSumDiffInclection=sign(areaSumDiff0.*areaSumDiff1);
areaSumDiffPos=sort(find(areaSumDiff1>0)); %pick the positive slope of the next point (indicating the valley point)
candidates0=sort(find(areaSumDiffInclection<0));
candidates=intersect(candidates0,areaSumDiffPos);
if ~isempty(candidates)
    areaBoundaryLength=candidates(end);
else
    areaBoundaryLength=5;
end
end