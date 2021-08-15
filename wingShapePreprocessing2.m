function [realCen, symAxis, symOrtho, boundingBox,  tipPts, refPts, wingParts, WingAxesSlopes, tformWingImg, shapeImg, allComImg, allinfo]=wingShapePreprocessing2(template,flag,mask,spp_json_directory,Result_directory,subFolderList,side,Sphingidae)
vdlist={'dorsal','ventral'};
%[specimenB,~]=bwboundaries(mask);
%sppEdgePt=specimenB{1};
%skeletonMask = bwmorph(mask,'skel',Inf);
% maskf=imfill(mask,'holes');
maskf=mask;
%maskf2=bwareafilt(imdilate(bwareafilt(imfill(imerode(mask,strel('disk',5)),'hole'),2),strel('disk',2)),1);%Used for reduced area

%Find the symetric axis
[symCentroid,symAxis,~]=findSymetricAxes(maskf);
disp('The symatric axis has been found.');
 %FInd the center based on regionprop function   
 regioncen0=regionprops(maskf,'Centroid','BoundingBox'); %The center of the bounding box
 regioncen= regioncen0.Centroid;
 boundingBox=regioncen0.BoundingBox;
 
 %Find the centroid of the eroded central region based on regionprop function   
 cenregionTrimPx=50;
 while 1
     cenregion=imerode(maskf, strel('disk', cenregionTrimPx)); %change to smaller value for small specimen
     if nnz(cenregion)>=5000
         break
     elseif cenregionTrimPx<=0
         break
     else
          cenregionTrimPx= cenregionTrimPx-10;
     end
 end
 cenregioncen0=regionprops(uint8(cenregion),'Centroid','BoundingBox'); %The center of the bounding box
 cenregioncen=cenregioncen0.Centroid;
 boundingBoxErosion=cenregioncen0.BoundingBox;

 %Calculated the difference between two centroids
 cenDiff=pdist([regioncen;cenregioncen],'euclidean');
 %cenDiffRatio=cenDiff/min(boundingBox(3:4));
 
 %Pick the most suitable one
 if cenDiff<50
    realCen= regioncen;
     disp('The centroid of the {entire mask} is used as the real cnetroid.');
 else
     realCen= cenregioncen;
      disp('The centroid of the {delated central region} is used as the real cnetroid.');
 end
 disp('The centroid has been determined.');
 
%Derive the coordination of corners
 ulCorner=boundingBox(1:2);
 urCorner=[boundingBox(1)+boundingBox(3),boundingBox(2)];
 llCorner=[boundingBox(1),boundingBox(2)+boundingBox(4)];
 lrCorner=[boundingBox(1)+boundingBox(3),boundingBox(2)+boundingBox(4)];
 allFrameCorners=[ ulCorner; urCorner; llCorner; lrCorner];
%%
%Prepare the symetric axes for ploting
%create symetric axes based on eigenvector
symOrtho=reshape(null(symAxis(:).'),1,[]);
dim_1=realCen+symAxis*size(maskf,1)/3;
dim_1plot=[realCen(1),dim_1(1);realCen(2),dim_1(2)];
dim_2=realCen+symOrtho*size(maskf,1)/9;
dim_2plot=[realCen(1),dim_2(1);realCen(2),dim_2(2)];

%figure,imshow(mask)
inspoutname=fullfile(Result_directory,'Shape_analysis',subFolderList{1},[template,'_',vdlist{side},flag,'_primary_key_pts.png']);
figinsp=figure('visible', 'off');
imshowpair(maskf,cenregion);
hold on;
plot(symCentroid(1),symCentroid(2),'y*');
plot(symCentroid(1),symCentroid(2),'bo');
%plot(symSkeletonCentroid(1),symSkeletonCentroid(2),'yo');
plot(regioncen(1),regioncen(2),'b*');
plot(cenregioncen(1), cenregioncen(2),'r+');
line(dim_1plot(1,:),dim_1plot(2,:), 'Color', 'r' ,'LineWidth', 1);
line(dim_2plot(1,:),dim_2plot(2,:), 'Color', 'b' ,'LineStyle','--','LineWidth', 1);
plot(allFrameCorners(:,1),allFrameCorners(:,2),'y+');
%rectangle('Position', boundingBox, 'EdgeColor','w','LineStyle','-.', 'LineWidth', 1);
%plot(realCen(1),realCen(2),'ro','LineWidth', 2);
hold off;
%saveas(figmask, maskoutname);
%print(figinsp,inspoutname,'-dpng','-r150'); %Use print to save high quality images
export_fig(figinsp,inspoutname, '-png','-r150');
close(figinsp);
disp('An image indcating the primary key points of specimen image has been saved.');
%%
%Show the mask in an image
maskoutname=fullfile(Result_directory,'Shape_analysis',subFolderList{2},[template,'_',vdlist{side},flag,'_mask_cen.png']);
figmask=figure('visible', 'off');
imshow(maskf);
hold on;
%plot(segPts(:,1),segPts(:,2),'r+');
%plot(forehindCornerR(:,1),forehindCornerR(:,2),'r*');
%plot(forehindCornerL(:,1),forehindCornerL(:,2),'r*');
plot(realCen(:,1),realCen(:,2),'bo','LineWidth', 1);
line(dim_1plot(1,:),dim_1plot(2,:), 'Color', 'b' ,'LineWidth', 1);
line(dim_2plot(1,:),dim_2plot(2,:), 'Color', 'b' ,'LineStyle','--','LineWidth', 1);
hold off;
%saveas(figmask, maskoutname);
%print(figmask,maskoutname,'-dpng','-r150'); %Use print to save high quality images
export_fig(figmask,maskoutname, '-png','-r150');
close(figmask);
disp('An image indcating the mask and centroid of specimen image has been saved.');
%%
%beltheight=boundingBox(4)*0.25;
%beltwidth=boundingBox(3)*0.15;
%Use erosion mask to prevent the interference of long tail
 if cenDiff>=50
    boundingBoxDV= boundingBoxErosion;
    ulCornerDV=boundingBoxErosion(1:2);
    lrCornerDV=[boundingBoxErosion(1)+boundingBoxErosion(3),boundingBoxErosion(2)+boundingBoxErosion(4)];   
 else
     boundingBoxDV= boundingBox;
     ulCornerDV=boundingBox(1:2);
     lrCornerDV=[boundingBox(1)+boundingBox(3),boundingBox(2)+boundingBox(4)];
 end
disp('########## Begin to find the corner between Fore and Hing wings. #########');
disp('Begin to find the corner between Left fore and hing wings.');
nStrongCornersList=[500,1000,2000,3000,4000]; %[500,1000,2000,4000];
nSectionList=[60:-5:20]; %number of elements should greater than 4  ([20:5:50];)

slopeSwitch='wingEdge';
[conjPt,forehindCorner,~]=findForeHindCorner(nStrongCornersList,nSectionList,maskf,realCen,symAxis,ulCornerDV,boundingBoxDV,slopeSwitch);
if length(forehindCorner(forehindCorner(:,1)>0))<5
    slopeSwitch='cenAxis';
    [conjPt,forehindCorner,~]=findForeHindCorner(nStrongCornersList,nSectionList,maskf,realCen,symAxis,ulCornerDV,boundingBoxDV,slopeSwitch);
end
forehindCornerL=conjPt;
disp('The corner between Left fore and hing wings has been found.');
disp('Begin to find the corner between Right fore and hing wings.');

slopeSwitch='wingEdge';
[conjPt,forehindCorner,~]=findForeHindCorner(nStrongCornersList,nSectionList,maskf,realCen,symAxis,lrCornerDV,boundingBoxDV,slopeSwitch);
if length(forehindCorner(forehindCorner(:,1)>0))<5
    slopeSwitch='cenAxis';
    [conjPt,forehindCorner,~]=findForeHindCorner(nStrongCornersList,nSectionList,maskf,realCen,symAxis,lrCornerDV,boundingBoxDV,slopeSwitch);
end
forehindCornerR=conjPt;
disp('The corner between Right fore and hing wings has been found.');
disp('########## Two corners between Fore and Hing wings are found. #########');
%%
%%
%Show the inspection image
maskoutname1=fullfile(Result_directory,'Shape_analysis',subFolderList{1},[template,'_',vdlist{side},flag,'_check_img.png']);
figmask1=figure('visible', 'off');
imshow(maskf);
hold on;
plot(forehindCornerR(:,1),forehindCornerR(:,2),'rx','LineWidth', 2);
plot(forehindCornerL(:,1),forehindCornerL(:,2),'rx','LineWidth', 2);
plot(realCen(:,1),realCen(:,2),'bo','LineWidth', 2);
line(dim_1plot(1,:),dim_1plot(2,:), 'Color', 'b' ,'LineWidth', 1);
line(dim_2plot(1,:),dim_2plot(2,:), 'Color', 'b' ,'LineStyle','--','LineWidth', 1);
hold off;
%saveas(figmask, maskoutname);
%print(figmask1,maskoutname1,'-dpng','-r150'); %Use print to save high quality images
export_fig(figmask1,maskoutname1, '-png','-r150');
close(figmask1);
disp('An image indcating the mask and fore hind wing corners has been saved.');
%%
%Find best beltWpar
disp('Start to find the optimal body width parameter.');
beltWpar=findbeltWpar(maskf,forehindCornerL,forehindCornerR,realCen,boundingBox);
disp('The body width parameter has been found.');
%%
disp('Start to crop upper and lower mask based on the fore-hing wing corners.');
upperSegMaskPts=[[0,0];[size(maskf,2),0];forehindCornerR;realCen;forehindCornerL];
upperSegMask = poly2mask(round(upperSegMaskPts(:,1)),round(upperSegMaskPts(:,2)),size(maskf,1),size(maskf,2));
upperMask=immultiply(maskf,upperSegMask);

lowerSegMaskPts=[[0,size(maskf,1)];[size(maskf,2),size(maskf,1)];forehindCornerR;realCen;forehindCornerL];
lowerSegMask = poly2mask(round(lowerSegMaskPts(:,1)),round(lowerSegMaskPts(:,2)),size(maskf,1),size(maskf,2));
lowerMask=immultiply(maskf,lowerSegMask);
disp('Upper and lower mask have been croped out based on the fore-hing wing corners.');

%%
disp('########## Start to find the corner between Wing and body. ##########');
disp('Begin to find the corner between Left-Fore wing and body.');
%Left fore wing
nStrongCornersList=[1000,1500,2000,3000,4000]; %[1000,1500,2000,4000];
%nSectionList=[20:4:40]; %number of elements should greater than 4 (lower section number is more conservative; finding shorter but better edge length)
nSectionList=reshape( [[20:4:50];[50:-4:20]],1,[]); %use a mix parameter lists low-high + high-low to searching for better results
[conjPt, conjCorners]=findBodyWingCorner2(nStrongCornersList,nSectionList,upperMask,realCen,symAxis,forehindCornerL,'LF',boundingBox,beltWpar,maskf);
conjCornerLF=conjPt;
disp('The corner between Left-Fore wing and body has been found.');

disp('Begin to find the corner between Right-Fore wing and body.');
%Right fore wing
[conjPt, conjCorners]=findBodyWingCorner2(nStrongCornersList,nSectionList,upperMask,realCen,symAxis,forehindCornerR,'RF',boundingBox,beltWpar,maskf);
conjCornerRF=conjPt;
disp('The corner between Right-Fore wing and body has been found.');

disp('Begin to find the corner between Left-Hind wing and body.');
%Left hind wing
nSectionList=[20:4:40]; %number of elements should greater than 4
[conjPt, conjCorners]=findBodyWingCorner2(nStrongCornersList,nSectionList,lowerMask,realCen,symAxis,forehindCornerL,'LH',boundingBox,beltWpar,maskf);
conjCornerLH=conjPt;
disp('The corner between Left-Hind wing and body has been found.');

disp('Begin to find the corner between Right-Hind wing and body.');
%Right hind wing
[conjPt, conjCorners]=findBodyWingCorner2(nStrongCornersList,nSectionList,lowerMask,realCen,symAxis,forehindCornerR,'RH',boundingBox,beltWpar,maskf);
conjCornerRH=conjPt;
disp('The corner between Right-Hind wing and body has been found.');
disp('########## All landmark corners are found. ##########');

%Integrate the landmark points: (L-ForeHind Corner, L-Fore-Body corner, R-Fore-Body corner, R-Fore-Hind Corner,R-Hind-Body corner, L-Hind-Body corner)
segPts=[forehindCornerL;conjCornerLF;conjCornerRF;forehindCornerR;conjCornerRH;conjCornerLH];


%%
%Find the slopes of wings' main axes: (LeftFore, LeftHind, RightFore, RightHind)
disp('########## Start to find the slope of long axes of wings. ##########');
[WingAxesSlopes,tipPts]=findSlopesOfWingAxes(maskf,symAxis,realCen,segPts,boundingBox);
disp('########## The slope of long axes of wings are determined. ##########');
%%%%%%%%%%%%%%%Special Region For Sphingidae%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For Sphingidae use this approach to derive the wingAxesSlopes (no hind wing rotation)
%Sphingidae=1; %1: Sphingidae; 0: not Sphingidae
if Sphingidae==1
    WingAxesSlopes(2)=0;
    WingAxesSlopes(4)=0;
    disp('The rotation of hind wings will not be corrected for Sphingidae.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%Special Region For Sphingidae%%%%%%%%%%%%%%%
%%
%Show characters in an image
charoutname=fullfile(Result_directory,'Shape_analysis',subFolderList{3},[template,'_',vdlist{side},flag,'_main_characters.png']);
figchar=figure('visible', 'off');
imshowpair(mask,cenregion);
hold on;
plot(symCentroid(1),symCentroid(2),'bo');
plot(symCentroid(1),symCentroid(2),'y+','LineWidth', 1);
%plot(symSkeletonCentroid(1),symSkeletonCentroid(2),'yo');
plot(regioncen(1),regioncen(2),'b*');
plot(cenregioncen(1), cenregioncen(2),'r+');
line(dim_1plot(1,:),dim_1plot(2,:), 'Color', 'r' ,'LineWidth', 1);
line(dim_2plot(1,:),dim_2plot(2,:), 'Color', 'b' ,'LineStyle','--','LineWidth', 1);
plot(segPts(:,1),segPts(:,2),'r+');
plot(forehindCornerR(:,1),forehindCornerR(:,2),'ro','LineWidth', 1);
plot(forehindCornerL(:,1),forehindCornerL(:,2),'ro','LineWidth', 1);
plot(realCen(:,1),realCen(:,2),'bo','LineWidth', 1);
plot(allFrameCorners(:,1),allFrameCorners(:,2),'y+');
plot(tipPts(:,1),tipPts(:,2),'yd','LineWidth', 1);
plot(tipPts(:,1),tipPts(:,2),'y+','LineWidth', 1);
hold off;
%saveas(figmask, maskoutname);
%print(figchar,charoutname,'-dpng','-r150'); %Use print to save high quality images
export_fig(figchar,charoutname, '-png','-r150');
close(figchar);
disp('An image indcating all key characters of specimen image has been saved.');
%%
%wingPartNameList={'Left Fore Wing','Left Hind Wing','Right Fore Wing','Right Hind Wing'};
%Segment the wings: (LeftFore, LeftHind, RightFore, RightHind)
disp('########## Start to chop four wing parts out of the image. ##########');
%find if there is a corresponding json file
json_ds = struct2dataset(dir(fullfile(spp_json_directory,[template,'_',vdlist{side},'*.json'])));
if ~isempty(json_ds)
    try
        disp('Find the corresponding json file');
        jsoninname=json_ds.name;
        json_data=loadjson(fullfile(spp_json_directory,jsoninname));
        [wingParts,refPts]=segWings_with_json3(mask,realCen,segPts,json_data);
        wingPartJudge=mask-wingParts{1}-wingParts{2}-wingParts{3}-wingParts{4};
        remainingRatio=nnz(wingPartJudge)/nnz(mask);
        if remainingRatio>0.3
            [wingParts,refPts]=segWings_with_json2(mask,realCen,segPts,json_data); %additional logic: pick one based on the position
        end
         disp('Complete with both mask and json file for segmentation');
    catch
        disp('Something goes wrong in the analysis with json file, return to use naive segmentation');
        [wingParts,refPts]=segWings(mask,realCen,segPts);
    end
else
    disp('Use only mask for segmentation');
    [wingParts,refPts]=segWings(mask,realCen,segPts);
end
disp('########## Four wing parts are cropped. ##########');
%figure,imshow(mask); hold on; plot(segPts(:,1),segPts(:,2),'r+');
%figure,imshow(mask); hold on; plot(refPts(:,1),refPts(:,2),'r+');
disp('########## Start to rotate each wing parts. ##########');
[tformWingImg, tformPtData]=wingRotation(wingParts,WingAxesSlopes,refPts);
disp('########## The rotation of all wings has been finished. ##########');

disp('########## Start to compose Fore wing and Hind wing. ##########');
%rotate and normalize the position of a pair of fore wing and hind wings
shapeImg=combineForeHindWings(tformWingImg,tformPtData);
disp('########## Fore wing and Hind wing have been composed. ##########');

disp('########## Start to integrate 4 wings into 1 panel. ##########');
[allComImg, allinfo]=integrate4WingsInto1Panel(shapeImg,tformWingImg,refPts);
disp('########## Integrated image has been generated. ##########');
%%
%Crop the blank region of image
[NewAllComImg,NewOrigion]=rmBlankRegion(allComImg,100);

%Show wing shape image
shapevisoutname=fullfile(Result_directory,'Shape_analysis',subFolderList{4},[template,'_',vdlist{side},flag,'_wings_shape.png']);
figshape=figure('visible', 'off');
imshow(NewAllComImg);
%saveas(figmask, maskoutname);
%print(figshape,shapevisoutname,'-dpng','-r150'); %Use print to save high quality images
export_fig(figshape,shapevisoutname, '-png','-r150');
close(figshape);
disp('An image showing adjusted specimen wings shape has been saved.');

shapevisoutname=fullfile(Result_directory,'Shape_analysis',subFolderList{4},[template,'_',vdlist{side},flag,'_wings_shape_n_pts.png']);
figshape=figure('visible', 'off');
imshow(NewAllComImg);
hold on;
plot(allinfo(:,1)-NewOrigion(1),allinfo(:,2)-NewOrigion(2),'yx','LineWidth', 2);
%line(allinfo(1:2,1)-NewOrigion(1),allinfo(1:2,2)-NewOrigion(2), 'Color', 'y' ,'LineWidth', 2);
%line(allinfo(3:4,1)-NewOrigion(1),allinfo(3:4,2)-NewOrigion(2), 'Color', 'y' ,'LineWidth', 2);
hold off;
%saveas(figmask, maskoutname);
%print(figshape,shapevisoutname,'-dpng','-r150'); %Use print to save high quality images
export_fig(figshape,shapevisoutname, '-png','-r150');
close(figshape);
disp('An image showing adjusted specimen wings shape and boundary Pts has been saved.');
end