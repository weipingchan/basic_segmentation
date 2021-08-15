function [bodyMask,bodyCharacters, antennaMask,antennaCharacters]=body_antenna_module2(mask,wingParts,refPts,tipPts,bodyTrimPx, antTrimPx, scalelen)
%%
% bodyTrimPx=5;
%Calculate morphological parameters
disp('Begin to calculate morphology of body and antenna');
%Body
wingMask=wingParts{1}+wingParts{2}+wingParts{3}+wingParts{4};
realCen=refPts(8,:); %rough antenna preparation
forehindCornerL=tipPts(1,:);
forehindCornerR=tipPts(2,:);
body0=bwareafilt(logical(mask-wingMask), 1);
body0Stat=regionprops(body0,'Centroid', 'MajorAxisLength','MinorAxisLength');
if antTrimPx > body0Stat.MinorAxisLength/3; antTrimPx= round(body0Stat.MinorAxisLength/3); end
majorObj=imdilate(bwareafilt(logical(imerode(mask,strel('disk', antTrimPx))), 1),strel('disk', antTrimPx)); %Used for reduced area
upperSegMaskPts=[[0,0];[size(mask,2),0];[size(mask,2),forehindCornerR(2)];forehindCornerR; refPts(3,:); refPts(2,:) ;forehindCornerL;[0,forehindCornerL(2)]];
upperSegMask = poly2mask(round(upperSegMaskPts(:,1)),round(upperSegMaskPts(:,2)),size(mask,1),size(mask,2));
antennaraw=mask-majorObj;
antenna=bwareafilt(logical(bwareaopen(immultiply(antennaraw,upperSegMask),200)),2);
bodyraw=mask-wingMask-antenna;
body=bwareafilt(logical(imdilate(imerode(bodyraw,strel('disk',bodyTrimPx)),strel('disk',bodyTrimPx))),1);
propBody = regionprops(body,'MinorAxisLength','MajorAxisLength');

%bodyInspect=body+(mask)*0.1;
if ~isempty(propBody)
    bodyLength=propBody.MajorAxisLength/scalelen;
    bodyWidth=propBody.MinorAxisLength/scalelen;
    bodyCharacters=[bodyLength,bodyWidth];
else
    bodyLength=-9999;
    bodyWidth=-9999;
    bodyCharacters=[bodyLength,bodyWidth];
end
bodyMask=body;

%Antenna %The unit for all antenna parameters is mm
minimalAntennaLength=60; %The length of antenna should longer than this length or will be neglected
 [antennaCharacters,antFig]=antennaMeasure8(mask,refPts,tipPts,wingMask,body,minimalAntennaLength, bodyraw, upperSegMask, scalelen);
 antennaMask=antFig>0.8;
 disp('Body and antenna are segmented');
end