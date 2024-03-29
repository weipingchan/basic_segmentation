function  tangentLineSlpoe=find_inner_tangent_line_slope2(partMask,symAxis,tarCorner,LeftRight, nSection,boundingBox)
%This funciton is only used for generating preliminary visualization
verVector=symAxis*max(size(partMask));
symOrtho=reshape(null(symAxis(:).'),1,[]);
horVector=symOrtho*max(size(partMask));

%Choose the corner of the bounding box
if LeftRight=='L'
    corn=[1, size(partMask,2)];
    verEdgeLine= [[1,1] ; corn];
    LeftRightForeHind='LF';
    boundingBoxCorner=[boundingBox(1), boundingBox(2)+boundingBox(4)];
else
    corn=flip(size(partMask));
    verEdgeLine=[[size(partMask,2),1] ; corn];
    LeftRightForeHind='RF';
    boundingBoxCorner=[boundingBox(1)+boundingBox(3), boundingBox(2)+boundingBox(4)];
end

horEdgeLine=[[1,size(partMask,1)] ; flip(size(partMask))];

verLine=[tarCorner-verVector; tarCorner+verVector];
horLine=[tarCorner-horVector; tarCorner+horVector];

[bintX,bintY] = polyxpoly(verLine(:,1),verLine(:,2),horEdgeLine(:,1),horEdgeLine(:,2));
[sintX,sintY] = polyxpoly(horLine(:,1),horLine(:,2),verEdgeLine(:,1),verEdgeLine(:,2));

intS=[sintX,sintY];
intB=[bintX,bintY];

ptGroup=[tarCorner; intS; corn; intB];
region=poly2mask(ptGroup(:,1),ptGroup(:,2),size(partMask,1),size(partMask,2));
wingMask=imdilate(imerode(immultiply(partMask,region),strel('disk',30)),strel('disk',29));

[upspecimenB,~]=bwboundaries(wingMask);
upperEdgePt=upspecimenB{1};
Corners=flip(upperEdgePt,2);

beltB=[boundingBoxCorner-horVector ; boundingBoxCorner+horVector];

usDL=tarCorner-boundingBoxCorner;
uSlopeL=usDL(2)/usDL(1);

Ls2cen=findSlopesForEvenAngle2(symAxis,uSlopeL,nSection,LeftRightForeHind);

%Calculate the intersection points
intersectAll=zeros(length(Ls2cen),2);
for slpn=1:length(Ls2cen)
    sVector=[1, Ls2cen(slpn)]*max(size(partMask));
    tmpSegPts=[tarCorner-sVector ; tarCorner+sVector];
    [intersectX,intersectY]= polyxpoly(tmpSegPts(:,1),tmpSegPts(:,2),beltB(:,1),beltB(:,2));
    intersectAll(slpn,:) = [intersectX,intersectY];
end

edgePtsCount=zeros(length(intersectAll),0);
for ccc=1:length(intersectAll)-1
    intersectPts=intersectAll(ccc:ccc+1,:);
    triRegion=[intersectPts;tarCorner];
    inPts =inpolygon(Corners(:,1),Corners(:,2),triRegion(:,1),triRegion(:,2));
    edgePtsCount(ccc)=sum(inPts);
end

edgePtsCountDiff=diff(edgePtsCount);
edgePtsCountDiff2=sign([edgePtsCountDiff,0]);
IdxLinear=findchangepts(edgePtsCount,'MaxNumChanges',3,'Statistic','linear');
IdxStd=findchangepts(edgePtsCount,'MaxNumChanges',3,'Statistic','std'); 
edgePtsCountDiff3=ones( [1,length(edgePtsCount)] );
edgePtsCountDiff3([IdxLinear,IdxStd]+1)=2; %point shows the disconnectivity to its neighbor
edgePtsCountMean=edgePtsCount;
edgePtsCountMean(edgePtsCountMean<mean(edgePtsCountMean))=0;

intersectLoc=find(edgePtsCountDiff2.*edgePtsCountDiff3.*sign(edgePtsCountMean)<=-2, 1 );
if isempty(intersectLoc)
    intersectLoc=find(edgePtsCountDiff3.*sign(edgePtsCountMean)>=2, 1 );
    if isempty(intersectLoc)
    intersectLoc=find(edgePtsCountDiff3.*sign(edgePtsCountMean)>=1, 1 );
    end
end

tangentLineSlpoe=Ls2cen(intersectLoc);

%Preserved for visualization when debugging
% vector= [1, tangentLineSlpoe]*size(partMask,1);
% vec=[tarCorner; tarCorner+sign(tangentLineSlpoe)*vector];
% figure,imshow(wingMask);hold on;
% plot(vec(:,1),vec(:,2),'r');
 
end