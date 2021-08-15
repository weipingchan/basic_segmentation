function [wingParts,refPts]=segWings_with_json2(mask,realCen,segPts,json_data)

%Retrive the line in original image size
[hind_fore_sep_l,resize_ratio]=get_original_line(mask,json_data.image,json_data.image.l_track);
[hind_fore_sep_r,~]=get_original_line(mask,json_data.image,json_data.image.r_track);

[l_row,l_col,~] = find(hind_fore_sep_l~=0);
coor_l=flip(sortrows([l_row,l_col],2,'descend'),2);
body_fh_pt_l=coor_l(1,:);

[r_row,r_col,~] = find(hind_fore_sep_r~=0);
coor_r=flip(sortrows([r_row,r_col],2,'ascend'),2);
body_fh_pt_r=coor_r(1,:);

% %Visualize json defined fore-hind wing seg line
% figure,imshow(mask);hold on;
% plot(coor_l(:,1),coor_l(:,2),'r');
% plot(coor_r(:,1),coor_r(:,2),'r');

%segPts=[forehindCornerL;conjCornerLF;conjCornerRF;forehindCornerR;conjCornerRH;conjCornerLH];
extendLength=max(size(mask));
LWingLineSeg=[segPts(2,:);segPts(6,:)];
%LWingVector=(LWingLineSeg(1,:)-LWingLineSeg(2,:))*max(size(partMask));
LWingLine=seg2line(LWingLineSeg,extendLength);

RWingLineSeg=[segPts(3,:);segPts(5,:)];
%RWingVector=(RWingLineSeg(1,:)-RWingLineSeg(2,:))*max(size(partMask));
RWingLine=seg2line(RWingLineSeg,extendLength);

LFHCenLineSeg=[segPts(1,:);realCen];
RFHCenLineSeg=[segPts(4,:);realCen];

LFcropLineSeg=segPts(1:2,:);
%LFcropVector=(LFcropLineSeg(1,:)-LFcropLineSeg(2,:))*max(size(partMask));
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

%Original setting
[intersectX,intersectY]= polyxpoly(LWingLine(:,1),LWingLine(:,2),LFHCenLineSeg(:,1),LFHCenLineSeg(:,2));
LWingRCornerPt=[intersectX,intersectY];
[intersectX,intersectY]= polyxpoly(RWingLine(:,1),RWingLine(:,2),RFHCenLineSeg(:,1),RFHCenLineSeg(:,2));
RWingLCornerPt=[intersectX,intersectY];

%Swap out with json manually define data
LWingRCornerPt=body_fh_pt_l;
RWingLCornerPt=body_fh_pt_r;

extendPx=10; %searching area for next step
disp(['Variable [extendPx]: ',num2str(extendPx)]);
pathL=find_wing_seg_path_line(mask, segPts(1,:), extendPx, 'L'); %Left
pathR=find_wing_seg_path_line(mask, segPts(4,:), extendPx, 'R'); %Right
pathBL=find_wing_seg_path_line(mask, segPts(6,:), extendPx, 'B'); %Bottom Left
pathBR=find_wing_seg_path_line(mask, segPts(5,:), extendPx, 'B'); %Bottom Right

%Show the cutting path
% figure,imshow(mask);hold on;
% plot(pathL(:,1),pathL(:,2),'r','lineWidth',2);
% plot(pathR(:,1),pathR(:,2),'r','lineWidth',2);
% plot(pathBL(:,1),pathBL(:,2),'r','lineWidth',2);
% plot(pathBR(:,1),pathBR(:,2),'r','lineWidth',2);

trimPx=3;
disp(['Variable [trimPx]: ',num2str(trimPx)]);
%Left Fore Wing
[intersectX,intersectY]= polyxpoly(LWingLine(:,1),LWingLine(:,2),Ubound(:,1),Ubound(:,2));
UboundPt=[intersectX,intersectY];
[intersectX,intersectY]= polyxpoly(LFcropLine(:,1),LFcropLine(:,2),Lbound(:,1),Lbound(:,2));
LboundPt=[intersectX,intersectY];
%LFcropPtSet=[1 1 ; UboundPt ; segPts(2,:) ; LWingRCornerPt ; segPts(1,:) ; LboundPt];
LFcropPtSet=[1 1 ; UboundPt ; segPts(2,:) ; LWingRCornerPt ;  coor_l ;segPts(1,:) ; pathL; LboundPt];
LFcropMask = roipoly(mask,LFcropPtSet(:,1),LFcropPtSet(:,2));
% LFwing=bwareafilt(logical(imdilate(imerode(immultiply(mask,LFcropMask),strel('disk',trimPx)),strel('disk',trimPx))),1);
LFwing=wing_selection_within_seg_with_json(mask, LFcropMask, trimPx, 'upper');
disp('Left Fore Wing  is cropped out.');

%Left Hind Wing
[intersectX,intersectY]= polyxpoly(LWingLine(:,1),LWingLine(:,2),Bbound(:,1),Bbound(:,2));
BboundPt=[intersectX,intersectY];
[intersectX,intersectY]= polyxpoly(LHcropLine(:,1),LHcropLine(:,2),Lbound(:,1),Lbound(:,2));
LboundPt=[intersectX,intersectY];
%LHcropPtSet=[1 size(mask,1) ; BboundPt ; segPts(6,:) ; LWingRCornerPt ; segPts(1,:) ; LboundPt];
LHcropPtSet=[1 size(mask,1) ; BboundPt ; flip(pathBL,1); segPts(6,:) ; LWingRCornerPt ;  coor_l ;segPts(1,:) ; pathL; LboundPt];
LHcropMask = roipoly(mask,LHcropPtSet(:,1),LHcropPtSet(:,2));
% LHwing=bwareafilt(logical(imdilate(imerode(immultiply(mask,LHcropMask),strel('disk',trimPx)),strel('disk',trimPx))),1);
LHwing=wing_selection_within_seg_with_json(mask, LHcropMask, trimPx, 'lower');
disp('Left Hind Wing  is cropped out.');

%Right Fore Wing
[intersectX,intersectY]= polyxpoly(RWingLine(:,1),RWingLine(:,2),Ubound(:,1),Ubound(:,2));
UboundPt=[intersectX,intersectY];
[intersectX,intersectY]= polyxpoly(RFcropLine(:,1),RFcropLine(:,2),Rbound(:,1),Rbound(:,2));
RboundPt=[intersectX,intersectY];
%RFcropPtSet=[size(mask,2) 1 ; UboundPt ; segPts(3,:) ; RWingLCornerPt ; segPts(4,:) ; RboundPt];
RFcropPtSet=[size(mask,2) 1 ; UboundPt ; segPts(3,:) ; RWingLCornerPt ; coor_r ; segPts(4,:) ; pathR; RboundPt];
RFcropMask = roipoly(mask,RFcropPtSet(:,1),RFcropPtSet(:,2));
% RFwing=bwareafilt(logical(imdilate(imerode(immultiply(mask,RFcropMask),strel('disk',trimPx)),strel('disk',trimPx))),1);
RFwing=wing_selection_within_seg_with_json(mask, RFcropMask, trimPx, 'upper');
disp('Right Fore Wing  is cropped out.'); 

%Right Hind Wing
[intersectX,intersectY]= polyxpoly(RWingLine(:,1),RWingLine(:,2),Bbound(:,1),Bbound(:,2));
BboundPt=[intersectX,intersectY];
[intersectX,intersectY]= polyxpoly(RHcropLine(:,1),RHcropLine(:,2),Rbound(:,1),Rbound(:,2));
RboundPt=[intersectX,intersectY];
%RHcropPtSet=[flip(size(mask)) ; BboundPt ; segPts(5,:) ; RWingLCornerPt ; segPts(4,:) ; RboundPt];
RHcropPtSet=[flip(size(mask)) ; BboundPt ; flip(pathBR,1); segPts(5,:) ; RWingLCornerPt ; coor_r ; segPts(4,:) ; pathR; RboundPt];
RHcropMask = roipoly(mask,RHcropPtSet(:,1),RHcropPtSet(:,2));
% RHwing=bwareafilt(logical(imdilate(imerode(immultiply(mask,RHcropMask),strel('disk',trimPx)),strel('disk',trimPx))),1);
RHwing=wing_selection_within_seg_with_json(mask, RHcropMask, trimPx, 'lower');
disp('Right Hind Wing  is cropped out.');



wingParts={smoothShape(LFwing), LHwing, smoothShape(RFwing), RHwing};
refPts=[segPts;LWingRCornerPt;realCen;RWingLCornerPt];

    function smoothResult=smoothShape(inMask)
        smoothResult=imdilate(bwareafilt(logical(imerode(inMask,strel('disk',5))),1),strel('disk',5));%Used for reduced area
    end
end