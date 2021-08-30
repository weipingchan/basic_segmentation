function conjCorners=find_body_wing_corner3(partMask,nStrongCorners,realCen,symAxis,tarCorner,LeftRightForeHind, nSection,boundingBox,beltWpar)
%Clean up the mask a little bit
trimPx0=5;
disp(['Variable [trimPx0]: ',num2str(trimPx0)]);
partMask0=bwareaopen(imdilate(imerode(partMask,strel('disk',trimPx0)),strel('disk',trimPx0)),300);

%Block the edge effect of the opposite side
if LeftRightForeHind=='LF'
    addMaskPts=[[size(partMask,2),0];[realCen+symAxis*max(size(partMask))];  [realCen-symAxis*max(size(partMask))]; [size(partMask,2),size(partMask,1)]];
elseif LeftRightForeHind=='RF'
    addMaskPts=[[0,0];[realCen+symAxis*max(size(partMask))];  [realCen-symAxis*max(size(partMask))]; [0,size(partMask,1)]];
elseif LeftRightForeHind=='LH'
    addMaskPts=[[size(partMask,2),0];[realCen+symAxis*max(size(partMask))];  [realCen-symAxis*max(size(partMask))]; [size(partMask,2),size(partMask,1)]];
else
    addMaskPts=[[0,0];[realCen+symAxis*max(size(partMask))];  [realCen-symAxis*max(size(partMask))]; [0,size(partMask,1)]];
end

addPartMask=poly2mask(round(addMaskPts(:,1)),round(addMaskPts(:,2)),size(partMask,1),size(partMask,2));
partMask1=imbinarize(partMask0+addPartMask,0.5);

%Find all edge points
[upspecimenB0,~]=bwboundaries(partMask1);
upperEdgePt0=upspecimenB0{1};
%%
%Detect all salinent points
corners = detectHarrisFeatures(partMask1);
if size(corners,1)<nStrongCorners*2/3
    Corners=corners.Location;   
elseif size(corners,1)<nStrongCorners
    Corners=corners.Location;    
else
    Corners=corners.selectStrongest(nStrongCorners).Location;
end
%%
%Calculate necessary vectors
verVector=symAxis*size(partMask,1);

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
        shortDIstPts=zeros(size(intersectPts,1),0);        
        for ppt=1:size(intersectPts,1)
            shortDIstPts(ppt)=pdist([intersectPts(ppt,:);tarCorner0],'euclidean');
        end
        sortedShortDIstPts= sort(shortDIstPts(shortDIstPts>emptyRegionLength));
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
%%
if LeftRightForeHind=='LF'
    conjCorners=leftconjCorners1(1,:);
elseif LeftRightForeHind=='RF'
    conjCorners=leftconjCorners1(1,:);
elseif LeftRightForeHind=='LH'
    conjCorners=findCloestPt(leftconjCorners1,realCen);
else
    conjCorners=findCloestPt(leftconjCorners1,realCen); %This is the original one
end

if ~isempty(conjCorners)
    disp('Find the key point.');
else
    disp('DID NOT find the key point.');
end
end