function pathList=find_wing_seg_path_line(mask, beginPt, extendPx, searchDirection)

% searchDirection='L'; %R B
% extendPx=10;
% beginPt=segPts(1,:);

pathList=beginPt;
while 1
reachRing=createCirclesMask(mask, pathList(end,:), extendPx+4)-createCirclesMask(mask, pathList(end,:), extendPx);
possibleMove=immultiply(reachRing, imcomplement(mask));
[moveX, moveY]=find(possibleMove);
possibleMovePts=[moveX, moveY];
if searchDirection=='L'
    possibleMovePts1=possibleMovePts(possibleMovePts(:,2)==min(possibleMovePts(:,2)),:);
    movePt=flip(possibleMovePts1(randsample(size(possibleMovePts1,1),1),:));
    pathList=[pathList; movePt];
    if movePt(1)==1
        break
    elseif size(pathList,1)>size(mask,2)/extendPx
        pathList=[];
        disp('CANNOT find a left seg-pathway');
        break
    end
elseif searchDirection=='R'
    possibleMovePts1=possibleMovePts(possibleMovePts(:,2)==max(possibleMovePts(:,2)),:);
    movePt=flip(possibleMovePts1(randsample(size(possibleMovePts1,1),1),:));
    pathList=[pathList; movePt];
    if movePt(1)==size(mask,2)
        break
    elseif size(pathList,1)>size(mask,2)/extendPx
        pathList=[];
        disp('CANNOT find a right seg-pathway');
        break
    end
elseif searchDirection=='B'
    possibleMovePts1=possibleMovePts(possibleMovePts(:,1)==max(possibleMovePts(:,1)),:);
    movePt=flip(possibleMovePts1(randsample(size(possibleMovePts1,1),1),:));
    pathList=[pathList; movePt];
    if movePt(2)==size(mask,1)
        break
    elseif size(pathList,1)>size(mask,1)/extendPx
        pathList=[];
        disp('CANNOT find a bottom seg-pathway');
        break
    end
end
end

end