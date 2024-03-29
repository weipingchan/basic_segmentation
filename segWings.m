function [wingParts,refPts]=segWings(mask,realCen,segPts)
%Use detected landmark on wing mask to segment wing pieces
%In this function, manually-defined fore-hindwing boundaries were not
%included; only solid lines are used
extendLength=max(size(mask));
LWingLineSeg=[segPts(2,:);segPts(6,:)];
LWingLine=seg2line(LWingLineSeg,extendLength);

RWingLineSeg=[segPts(3,:);segPts(5,:)];
RWingLine=seg2line(RWingLineSeg,extendLength);

LFHCenLineSeg=[segPts(1,:);realCen];
RFHCenLineSeg=[segPts(4,:);realCen];

LFcropLineSeg=segPts(1:2,:);
LFcropLine=seg2line(LFcropLineSeg,extendLength);

LHcropLineSeg=[segPts(6,:);segPts(1,:)];
LHcropLine=seg2line(LHcropLineSeg,extendLength);

RFcropLineSeg=segPts(3:4,:);
RFcropLine=seg2line(RFcropLineSeg,extendLength);

RHcropLineSeg=segPts(4:5,:);
RHcropLine=seg2line(RHcropLineSeg,extendLength);

Lbound=[1 1-1000*size(mask,1); 1 1000*size(mask,1)];
Rbound=[size(mask,2) 1-1000*size(mask,1); size(mask,2) 1000*size(mask,1)];
Ubound=[1-3*size(mask,2) 1; 3*size(mask,2) 1];
Bbound=[1-3*size(mask,2) size(mask,1); 3*size(mask,2) size(mask,1)];

[intersectX,intersectY]= polyxpoly(LWingLine(:,1),LWingLine(:,2),LFHCenLineSeg(:,1),LFHCenLineSeg(:,2));
LWingRCornerPt=[intersectX,intersectY];
[intersectX,intersectY]= polyxpoly(RWingLine(:,1),RWingLine(:,2),RFHCenLineSeg(:,1),RFHCenLineSeg(:,2));
RWingLCornerPt=[intersectX,intersectY];


trimPx=3;
disp(['Variable [trimPx]: ',num2str(trimPx)]);
%Left forewing
[intersectX,intersectY]= polyxpoly(LWingLine(:,1),LWingLine(:,2),Ubound(:,1),Ubound(:,2));
UboundPt=[intersectX,intersectY];
[intersectX,intersectY]= polyxpoly(LFcropLine(:,1),LFcropLine(:,2),Lbound(:,1),Lbound(:,2));
LboundPt=[intersectX,intersectY];
LFcropPtSet=[1 1 ; UboundPt ; segPts(2,:) ; LWingRCornerPt ; segPts(1,:) ; LboundPt];
LFcropMask = roipoly(mask,LFcropPtSet(:,1),LFcropPtSet(:,2));
LFwing=bwareafilt(logical(imdilate(imerode(immultiply(mask,LFcropMask),strel('disk',trimPx)),strel('disk',trimPx))),1);
disp('Left forewing is cropped out.');

%Left hindwing
[intersectX,intersectY]= polyxpoly(LWingLine(:,1),LWingLine(:,2),Bbound(:,1),Bbound(:,2));
BboundPt=[intersectX,intersectY];
[intersectX,intersectY]= polyxpoly(LHcropLine(:,1),LHcropLine(:,2),Lbound(:,1),Lbound(:,2));
LboundPt=[intersectX,intersectY];
LHcropPtSet=[1 size(mask,1) ; BboundPt ; segPts(6,:) ; LWingRCornerPt ; segPts(1,:) ; LboundPt];
LHcropMask = roipoly(mask,LHcropPtSet(:,1),LHcropPtSet(:,2));
LHwing=bwareafilt(logical(imdilate(imerode(immultiply(mask,LHcropMask),strel('disk',trimPx)),strel('disk',trimPx))),1);
disp('Left hindwing is cropped out.');

%Right forewing
[intersectX,intersectY]= polyxpoly(RWingLine(:,1),RWingLine(:,2),Ubound(:,1),Ubound(:,2));
UboundPt=[intersectX,intersectY];
[intersectX,intersectY]= polyxpoly(RFcropLine(:,1),RFcropLine(:,2),Rbound(:,1),Rbound(:,2));
RboundPt=[intersectX,intersectY];
RFcropPtSet=[size(mask,2) 1 ; UboundPt ; segPts(3,:) ; RWingLCornerPt ; segPts(4,:) ; RboundPt];
RFcropMask = roipoly(mask,RFcropPtSet(:,1),RFcropPtSet(:,2));
RFwing=bwareafilt(logical(imdilate(imerode(immultiply(mask,RFcropMask),strel('disk',trimPx)),strel('disk',trimPx))),1);
disp('Right forewing is cropped out.'); 

%Right hindwing
[intersectX,intersectY]= polyxpoly(RWingLine(:,1),RWingLine(:,2),Bbound(:,1),Bbound(:,2));
BboundPt=[intersectX,intersectY];
[intersectX,intersectY]= polyxpoly(RHcropLine(:,1),RHcropLine(:,2),Rbound(:,1),Rbound(:,2));
RboundPt=[intersectX,intersectY];
RHcropPtSet=[flip(size(mask)) ; BboundPt ; segPts(5,:) ; RWingLCornerPt ; segPts(4,:) ; RboundPt];
RHcropMask = roipoly(mask,RHcropPtSet(:,1),RHcropPtSet(:,2));
RHwing=bwareafilt(logical(imdilate(imerode(immultiply(mask,RHcropMask),strel('disk',trimPx)),strel('disk',trimPx))),1);
disp('Right hindwing is cropped out.');



wingParts={smoothShape(LFwing), LHwing, smoothShape(RFwing), RHwing};
refPts=[segPts;LWingRCornerPt;realCen;RWingLCornerPt];

    function smoothResult=smoothShape(inMask)
        smoothResult=imdilate(bwareafilt(logical(imerode(inMask,strel('disk',5))),1),strel('disk',5));%Used for reduced area
    end
end