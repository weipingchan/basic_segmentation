function [WingAxesSlopes,tipPts]=findSlopesOfWingAxes(mask,symAxis,realCen,segPts,boundingBox)

lowerSegMaskPts=[[0,size(mask,1)];[size(mask,2),size(mask,1)];segPts(4,:);realCen;segPts(1,:)];
lowerSegMask = roipoly(mask,round(lowerSegMaskPts(:,1)),round(lowerSegMaskPts(:,2)));
lowerMask=immultiply(mask,lowerSegMask);

skeletonMask = bwmorph(mask,'skel',Inf);
skeletonEndPtsMap = bwmorph(skeletonMask,'endpoint',Inf);
[skeletonEndPtX,skeletonEndPtY] = find(skeletonEndPtsMap==1);
skeletonEndPts=[skeletonEndPtY,skeletonEndPtX];
disp('The skeleton Image has been created.');
%figure,imshow(mask); hold on;
%plot(skeletonEndPts(:,1),skeletonEndPts(:,2),'r*');

%To prevent the result from round shape wing may be biased by thelateral
%edge, the front end of bounding box are lifted alittle bit.
liftedLength=100;
disp(['Variable [liftedLength]: ',num2str(liftedLength)]);
%Left end
tipL=findCloestPt(skeletonEndPts,boundingBox(1:2)-[0 liftedLength]);
%Right end
tipR=findCloestPt(skeletonEndPts,[boundingBox(1)+boundingBox(3),boundingBox(2)-liftedLength]);
%figure,imshow(mask);hold on;plot(tipL(1),tipL(2),'r*');plot(tipR(1),tipR(2),'r*');plot(segPts(2,1),segPts(2,2),'r*');
disp('Two tips of fore wing front ends are found.');

VectorL=tipL-segPts(2,:);
foreWIngLongAxisSlopeL=VectorL(2)/VectorL(1);
VectorR=tipR-segPts(3,:);
foreWIngLongAxisSlopeR=VectorR(2)/VectorR(1);

nSection=10;
partMask=lowerMask;
%Left
%tarCorner=segPts(end,:);
%LeftRight='L';
disp('Start to find the inner tangent line of Left hind wing.');
try
    tangentLineSlpoeL=find_inner_tangent_line_slope2(partMask,symAxis,segPts(end,:),'L', nSection,boundingBox);
    if isempty(tangentLineSlpoeL) tangentLineSlpoeL=symAxis(2)/symAxis(1);, end; 
catch
    tangentLineSlpoeL=symAxis(2)/symAxis(1);
end
disp('The slope of inner tangent line of Left hind wing is determined.');
%Right
disp('Start to find the inner tangent line of Rt hind wing.');
try
    tangentLineSlpoeR=find_inner_tangent_line_slope2(partMask,symAxis,segPts(end-1,:),'R', nSection,boundingBox);
    if isempty(tangentLineSlpoeR) tangentLineSlpoeR=symAxis(2)/symAxis(1);, end;
catch
    tangentLineSlpoeR=symAxis(2)/symAxis(1);
end
disp('The slope of inner tangent line of Right hind wing is determined.');

tipPts=[tipL; tipR];
WingAxesSlopes=[foreWIngLongAxisSlopeL, tangentLineSlpoeL, foreWIngLongAxisSlopeR, tangentLineSlpoeR];
end