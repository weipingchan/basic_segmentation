function conjCorners=find_body_wing_corner4(partMask,nStrongCorners,realCen,symAxis,tarCorner,LeftRightForeHind, nSection,boundingBox,beltWpar,maskf)

%Clean up the mask a little bit
trimPx0=5;
disp(['Variable [trimPx0]: ',num2str(trimPx0)]);
partMask0=bwareaopen(imdilate(imerode(partMask,strel('disk',trimPx0)),strel('disk',trimPx0)),300);

%Find all edge points for partMask
[upspecimenB00,~]=bwboundaries(partMask0);
upperEdgePt00=upspecimenB00{1};
%%
%Calculate necessary vectors
 if symAxis(2)<0 symAxis=-symAxis;, end;
verVector=symAxis*max(size(partMask));

symOrtho=reshape(null(symAxis(:).'),1,[]);
 if symOrtho(1)<0 symOrtho=-symOrtho;, end;
horVector=symOrtho*max(size(partMask));

%%
%Block the edge effect of the opposite side
ocut=[realCen+verVector;  realCen-verVector];
[ointsX,ointsY]= polyxpoly(ocut(:,1),ocut(:,2),upperEdgePt00(:,2),upperEdgePt00(:,1));
cutPts0=[ointsX,ointsY];

if LeftRightForeHind=='LF'
    addMaskPts=[[size(partMask,2),0];[realCen-verVector];  [realCen+verVector]; [size(partMask,2),size(partMask,1)]];
    cutPts1=cutPts0(cutPts0(:,2)<realCen(2)-10,:);
    cutDist = sqrt(sum(bsxfun(@minus,cutPts1,realCen).^2,2));
    cutPt=cutPts1(cutDist==min(cutDist),:);
    deMaskPts=[[size(partMask,2),0] ; [realCen-verVector]; cutPt ; [cutPt+horVector]];
elseif LeftRightForeHind=='RF'
    addMaskPts=[[0,0];[realCen-verVector];  [realCen+verVector]; [0,size(partMask,1)]];
    cutPts1=cutPts0(cutPts0(:,2)<realCen(2)-10,:);
    cutDist = sqrt(sum(bsxfun(@minus,cutPts1,realCen).^2,2));
    cutPt=cutPts1(cutDist==min(cutDist),:);
    deMaskPts=[[0,0];[realCen-verVector]; cutPt ; [cutPt-horVector]];
elseif LeftRightForeHind=='LH'
    addMaskPts=[[size(partMask,2),0];[realCen-verVector];  [realCen+verVector]; [size(partMask,2),size(partMask,1)]];
    cutPts1=cutPts0(cutPts0(:,2)>realCen(2)+10,:);
    cutDist = sqrt(sum(bsxfun(@minus,cutPts1,realCen).^2,2));
    cutPt=cutPts1(cutDist==min(cutDist),:);
    deMaskPts=[[size(partMask,2),size(partMask,1)] ; [realCen+verVector]; cutPt ; [cutPt+horVector]];
else
    addMaskPts=[[0,0];[realCen-verVector];  [realCen+verVector]; [0,size(partMask,1)]];
    cutPts1=cutPts0(cutPts0(:,2)>realCen(2)+10,:);
    cutDist = sqrt(sum(bsxfun(@minus,cutPts1,realCen).^2,2));
    cutPt=cutPts1(cutDist==min(cutDist),:);
    deMaskPts=[[0,size(partMask,1)];[realCen+verVector]; cutPt ; [cutPt-horVector]];
end

addPartMask=poly2mask(round(addMaskPts(:,1)),round(addMaskPts(:,2)),size(partMask,1),size(partMask,2));
dePartMask=poly2mask(round(deMaskPts(:,1)),round(deMaskPts(:,2)),size(partMask,1),size(partMask,2));
partMask1=imbinarize(partMask0+addPartMask,0.5);
partMask2=bwareaopen(imdilate(imerode(imbinarize(partMask1-dePartMask,0.5),strel('disk',trimPx0)),strel('disk',trimPx0)),300);
% figure,imshow(partMask2);

%Find all edge points for partMask
[upspecimenB0,~]=bwboundaries(partMask2);
upperEdgePt0=upspecimenB0{1};

%Find all edge points
[B0,~]=bwboundaries(maskf);
edgePt=B0{1};
%%
%Detect all salinent points
corners = detectHarrisFeatures(partMask2);
if size(corners,1)<nStrongCorners*2/3 %Legacy line
    Corners=corners.Location;   
elseif size(corners,1)<nStrongCorners
    Corners=corners.Location;    
else
    Corners=corners.selectStrongest(nStrongCorners).Location;
end

%%
%Derive the belt information; the belt should include the left and right fore-hindwing gaps for automatic detection
beltwidth=boundingBox(3)*beltWpar;

beltL=[realCen+[+round(beltwidth/2) 0]-verVector ; realCen+[+round(beltwidth/2) 0]+verVector];
beltR=[realCen+[-round(beltwidth/2) 0]-verVector ; realCen+[-round(beltwidth/2) 0]+verVector];
verBeltRegion=[beltL;flip(beltR,1)];
%%
emptyRegionLength=50;

%%
%Calculate the slope of tarCorner and centroid
usDL=tarCorner-realCen;
uSlopeL=usDL(2)/usDL(1);
%%
%move the reference point into the maks in order to prevent the interference of the zig-zag boundary
if LeftRightForeHind=='LF'
    tarCorner0=tarCorner+10*[1 -1];
elseif LeftRightForeHind=='RF'
    tarCorner0=tarCorner+10*[-1 -1];
elseif LeftRightForeHind=='LH'
    tarCorner0=tarCorner+10*[1 1];
else
    tarCorner0=tarCorner+10*[-1 1];
end

%%
%Derive slopes of all evenly spaced angle
Ls2cen=findSlopesForEvenAngle2(symAxis,uSlopeL,nSection,LeftRightForeHind);
%%
%Calculate the intersection points
intersectAll=cell(length(Ls2cen),0);
for slpn=1:length(Ls2cen)
    sVector=[1, Ls2cen(slpn)]*max(size(partMask));
    tmpSegPts=[tarCorner0-sVector ; tarCorner0+sVector]; %Move the point a little bit to prevent the integer coordination issue
    [intersectX,intersectY]= polyxpoly(tmpSegPts(:,1),tmpSegPts(:,2),upperEdgePt0(:,2),upperEdgePt0(:,1));
    intersectAll{slpn} = [intersectX,intersectY];
end
%%
%Calculate 3 indices for determine the targeted point
%1. Number of segments
%2. If there is candidate points in the belt region
%3. The abrupt change in the dinstance to the nearest point (wing edge)
intersectSegCount=zeros(length(intersectAll),0);
intersectDistPtsCount=zeros(length(intersectAll),0);
intersectShortestDist=zeros(length(intersectAll),0);
intersectDistPts=[];
for ccc=1:length(intersectAll)
    intersectSegCount(ccc)=length(intersectAll{ccc})/2;
    intersectPts=intersectAll{ccc};
    inPts0 =inpolygon(intersectPts(:,1),intersectPts(:,2),verBeltRegion(:,1),verBeltRegion(:,2));
    inPts=intersectPts(inPts0,:);
    intersectDistPtsCount(ccc)=size(inPts,1);
    cIdx=ccc+zeros(size(inPts,1), 1);
    intersectDistPts=cat(1,intersectDistPts,cat(2,cIdx,inPts));            
    
    %Calculate the distance of the nearest point
    if size(intersectPts,1)>0
       %emptyPts=inpolygon(intersectPts(:,1),intersectPts(:,2),emptyRegion(:,1),emptyRegion(:,2));
        %intersectPts2=intersectPts(~emptyPts,:);
        shortDIstPts=zeros(size(intersectPts,1),0);        
        for ppt=1:size(intersectPts,1)
            shortDIstPts(ppt)=pdist([intersectPts(ppt,:);tarCorner0],'euclidean');
        end
        sortedShortDIstPts= sort(shortDIstPts(shortDIstPts>emptyRegionLength));
        %sortedShortDIstPts = sort(shortDIstPts);
        %intersectShortestDist(ccc) = sortedShortDIstPts(2);
        intersectShortestDist(ccc) = sortedShortDIstPts(1); %The shortest distance
    else
        intersectShortestDist(ccc) = 0;
    end    
end

%If there is candidate points in the belt region
intersectDistPtsCount2=intersectDistPtsCount;
intersectDistPtsCount2(intersectDistPtsCount2>0)=1; %having value in belt -> 1

% Detect the changing point in the dinstance to the nearest point (wing edge)
intersectShortestDistDiff=diff(intersectShortestDist);
intersectShortestDistDiff2=sign([intersectShortestDistDiff,0]);
IdxLinear=findchangepts(intersectShortestDist,'MaxNumChanges',6,'Statistic','linear'); %THE NUMBER OF CHANGING PT here is sensitve to damaged wings. 6 is enough in most cases.
IdxStd=findchangepts(intersectShortestDist,'MaxNumChanges',6,'Statistic','std'); 
intersectShortestDistDiff3=ones( [1,length(intersectShortestDist)] );
intersectShortestDistDiff3([IdxLinear,IdxStd]-1)=2; %points shows the disconectivity to its neighbor
intersectShortestDistDiff4=intersectShortestDistDiff3;
intersectShortestDistDiff4(intersectShortestDistDiff2.*intersectShortestDistDiff4<0)=0;

% Detect the changing point of Number of segments
intersectSegCountDiff=diff(intersectSegCount);
intersectSegCountDiff2=[intersectSegCountDiff,0];
intersectSegCountDiff2(intersectSegCountDiff2<0)=-1; %The decrease of  segments represents the end of the front end edge

%Use all three indicies to determine the target segment line
intersectLocList=find(intersectDistPtsCount2.*intersectSegCountDiff2.*intersectShortestDistDiff4<=-2);

blockPts=cell(0,1);
ceid=1;
for locn=1:length(intersectLocList)
    intersectLoc=intersectLocList(locn);
    %Calculate the vector for final local searching
    block0=[tarCorner0-[1, Ls2cen(intersectLoc)]*max(size(partMask)) ; tarCorner0+[1, Ls2cen(intersectLoc)]*max(size(partMask))];
    block1=[tarCorner0-[1, Ls2cen(intersectLoc+1)]*max(size(partMask)) ; tarCorner0+[1, Ls2cen(intersectLoc+1)]*max(size(partMask))];

    %Calculate the four corners for final local searching
    try
        fineBox=zeros(4,2);
        [fineBox(1,1),fineBox(1,2)]=polyxpoly(block0(:,1),block0(:,2),beltL(:,1),beltL(:,2));
        [fineBox(2,1),fineBox(2,2)]=polyxpoly(block0(:,1),block0(:,2),beltR(:,1),beltR(:,2));
        [fineBox(3,1),fineBox(3,2)]=polyxpoly(block1(:,1),block1(:,2),beltR(:,1),beltR(:,2));
        [fineBox(4,1),fineBox(4,2)]=polyxpoly(block1(:,1),block1(:,2),beltL(:,1),beltL(:,2));
    catch %if segments are not long enough to have intersection, elongate the segments and try again
        beltL=[realCen+[+round(beltwidth/2) 0]-verVector*10^5 ; realCen+[+round(beltwidth/2) 0]+verVector*10^5];
        beltR=[realCen+[-round(beltwidth/2) 0]-verVector*10^5 ; realCen+[-round(beltwidth/2) 0]+verVector*10^5];

        block0=[tarCorner0-[1, Ls2cen(intersectLoc)]*max(size(partMask))*10^5 ; tarCorner0+[1, Ls2cen(intersectLoc)]*max(size(partMask))*10^5];
        block1=[tarCorner0-[1, Ls2cen(intersectLoc+1)]*max(size(partMask))*10^5 ; tarCorner0+[1, Ls2cen(intersectLoc+1)]*max(size(partMask))*10^5];

        fineBox=zeros(4,2);
        [fineBox(1,1),fineBox(1,2)]=polyxpoly(block0(:,1),block0(:,2),beltL(:,1),beltL(:,2));
        [fineBox(2,1),fineBox(2,2)]=polyxpoly(block0(:,1),block0(:,2),beltR(:,1),beltR(:,2));
        [fineBox(3,1),fineBox(3,2)]=polyxpoly(block1(:,1),block1(:,2),beltR(:,1),beltR(:,2));
        [fineBox(4,1),fineBox(4,2)]=polyxpoly(block1(:,1),block1(:,2),beltL(:,1),beltL(:,2));
    end

     %Pick up all points in the final searching box
    inCorner = inpolygon(Corners(:,1),Corners(:,2),fineBox(:,1),fineBox(:,2));
    leftconjCorners0=Corners(inCorner,:);
    
    if ~isempty(leftconjCorners0)
        blockPts{ceid}=leftconjCorners0;
        ceid=ceid+1;
    end
end

%This script is preserved for visualization when debugging
%     figure,imshow(partMask);hold on;
%     plot(block0(:,1),block0(:,2),'r')
%     plot(block1(:,1),block0(:,1),'r')
%     plot(fineBox(:,1),fineBox(:,2),'y')
%     plot(beltL(:,1),beltL(:,2),'r')
%     plot(beltR(:,1),beltR(:,2),'r')
%     plot(Corners(:,1),Corners(:,2),'g*')

%Pick the group that closest to the realCen
canPtlist=[];
for gid=1:length(blockPts)
    avgPt=mean(blockPts{gid});
    canPtlist=[canPtlist; avgPt];
end

groupid=find(ismember(canPtlist,findCloestPt(canPtlist,realCen),'rows'),1);

leftconjCorners1=blockPts{groupid};

%Remove those points within the mask
[inMask,~] = inpolygon(leftconjCorners1(:,1),leftconjCorners1(:,2),upperEdgePt00(:,2),upperEdgePt00(:,1));
leftconjCorners1(inMask,:)=[];
%%
%Calculate the distance to symOrtho at realCen
d2orthoVec=point_to_line_distance(leftconjCorners1, realCen+symOrtho, realCen-symOrtho);

%Calculate the distance to the opposite side of a specimen mask
if LeftRightForeHind=='LF'
    if symAxis(2)<0 symAxis=-symAxis;, end;
elseif LeftRightForeHind=='RF'
    if symAxis(2)<0 symAxis=-symAxis;, end;
elseif LeftRightForeHind=='LH'
    if symAxis(2)>0 symAxis=-symAxis;, end;
else
    if symAxis(2)>0 symAxis=-symAxis;, end;
end

symCutLineEndPt=leftconjCorners1+symAxis*max(size(partMask));

internalDistList=zeros(0,size(leftconjCorners1,1));
for segID=1:size(leftconjCorners1,1)
    startPt=leftconjCorners1(segID,:);
    endPt=symCutLineEndPt(segID,:);
    cutPts=[startPt+symAxis*2; endPt];
    
    [intsX,intsY]= polyxpoly(cutPts(:,1),cutPts(:,2),edgePt(:,2),edgePt(:,1));
    candEndPts=[intsX,intsY];
    travelDist = sqrt(sum(bsxfun(@minus,candEndPts,startPt).^2,2));
    internalDist=min(travelDist(travelDist>=pdist([tarCorner0;realCen])/2));
    if isempty(internalDist) internalDist=max(travelDist);, end;
    internalDistList(segID)=internalDist;
end

%Pick the one having shortest distance to the symOrtho at realCen and the shortest distance to the opposite side of a specimen mask
d2orthoVecIdx=d2orthoVec==min(d2orthoVec); %strict rule

internalDistList2=sort(internalDistList);
if length(internalDistList2)>=3 %Give top 3 best candidates
    internalIdx=ismember(internalDistList,internalDistList2(1:3));
else
    internalIdx=ismember(internalDistList,internalDistList2(1:end));
end

tarIdx=d2orthoVecIdx.*internalIdx';

if ~isempty(find(tarIdx))
    conjCorners0=leftconjCorners1(tarIdx,:);
else %If no consensus in two index, then use the first rule
    conjCorners0=leftconjCorners1(d2orthoVecIdx,:);
end

conjCorners=findCloestPt([edgePt(:,2),edgePt(:,1)],conjCorners0); %Find the corresponding Pt on the real edge

if ~isempty(conjCorners)
    disp('Find the key point.');
else
    disp('DID NOT find the key point.');
end
end