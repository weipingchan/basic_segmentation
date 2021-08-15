function [conjPt, conjCorners]=findBodyWingCorner2(nStrongCornersList,nSectionList,partMask,realCen,symAxis,tarCorner,LeftRightForeHind,boundingBox,beltWpar,maskf)

    [conjPt01, conjCorners01]=findBodyWingCorner01(nStrongCornersList,nSectionList,partMask,realCen,symAxis,tarCorner,LeftRightForeHind,boundingBox,beltWpar,maskf);

    [conjPt02, conjCorners02]=findBodyWingCorner02(nStrongCornersList,nSectionList,partMask,realCen,symAxis,tarCorner,LeftRightForeHind,boundingBox,beltWpar,maskf);
    %The function embeded here uses additional logic to pick the point
    %closer to the horizontal line at cnetroid; this logic is also applied
    %below
    
    conjPtAll=[conjPt01; conjPt02];
    conjCornersAll={conjCorners01,conjCorners02};
    
    symOrtho=reshape(null(symAxis(:).'),1,[]);
    
    cenOrtho=[realCen+symOrtho; realCen-symOrtho];
    
    %Calculate the distance to symOrtho at realCen
    d2orthoVec=point_to_line_distance(conjPtAll, cenOrtho(1,:),cenOrtho(2,:));
    optIdx=find(d2orthoVec==min(d2orthoVec));
    
    conjPt=conjPtAll(optIdx,:);
    conjCorners=conjCornersAll{optIdx};
    
%         figure,imshow(partMask);hold on;
%         plot(cenOrtho(:,1),cenOrtho(:,2),'r');
%         plot(conjPt(:,1),conjPt(:,2),'rO');
%         plot(realCen(:,1),realCen(:,2),'b*');
end