function [WingAxesSlopes,tipPts]=findSlopesOfWingAxes(mask,symAxis,realCen,segPts,boundingBox)
%This function is used for the preliminary visualization only
lowerSegMaskPts=[[0,size(mask,1)];[size(mask,2),size(mask,1)];segPts(4,:);realCen;segPts(1,:)];
lowerSegMask = roipoly(mask,round(lowerSegMaskPts(:,1)),round(lowerSegMaskPts(:,2)));
lowerMask=immultiply(mask,lowerSegMask);

skeletonMask = bwmorph(mask,'skel',Inf);
skeletonEndPtsMap = bwmorph(skeletonMask,'endpoint',Inf);
[skeletonEndPtX,skeletonEndPtY] = find(skeletonEndPtsMap==1);
skeletonEndPts=[skeletonEndPtY,skeletonEndPtX];
disp('The skeleton Image has been created.');
%The script is preserved for visualization when debugging
%figure,imshow(mask); hold on;
%plot(skeletonEndPts(:,1),skeletonEndPts(:,2),'r*');

%To prevent the result where round-shaped wing may be biased by the lateral
%edge, the front ends of bounding box are lifted a little bit.
liftedLength=100;
disp(['Variable [liftedLength]: ',num2str(liftedLength)]);
%Left end
tipL=findCloestPt(skeletonEndPts,boundingBox(1:2)-[0 liftedLength]);
%Right end
tipR=findCloestPt(skeletonEndPts,[boundingBox(1)+boundingBox(3),boundingBox(2)-liftedLength]);
disp('Two tips of forewing front ends are found.');

VectorL=tipL-segPts(2,:);
foreWIngLongAxisSlopeL=VectorL(2)/VectorL(1);
VectorR=tipR-segPts(3,:);
foreWIngLongAxisSlopeR=VectorR(2)/VectorR(1);

nSection=10;
partMask=lowerMask;
disp('Start to find the inner tangent line of the left hindwing.');
try
    tangentLineSlpoeL=find_inner_tangent_line_slope2(partMask,symAxis,segPts(end,:),'L', nSection,boundingBox);
    if isempty(tangentLineSlpoeL) tangentLineSlpoeL=symAxis(2)/symAxis(1);, end; 
catch
    tangentLineSlpoeL=symAxis(2)/symAxis(1);
end
disp('The slope of the inner tangent line of the left hindwing is determined.');
%Right
disp('Start to find the inner tangent line of the right hindwing.');
try
    tangentLineSlpoeR=find_inner_tangent_line_slope2(partMask,symAxis,segPts(end-1,:),'R', nSection,boundingBox);
    if isempty(tangentLineSlpoeR) tangentLineSlpoeR=symAxis(2)/symAxis(1);, end;
catch
    tangentLineSlpoeR=symAxis(2)/symAxis(1);
end
disp('The slope of the inner tangent line of the right hindwing is determined.');

tipPts=[tipL; tipR];
WingAxesSlopes=[foreWIngLongAxisSlopeL, tangentLineSlpoeL, foreWIngLongAxisSlopeR, tangentLineSlpoeR];
end