function [realCen, symAxis, symOrtho, boundingBox,  tipPts, refPts, wingParts, WingAxesSlopes, tformWingImg, shapeImg, allComImg, allinfo]=wingShapePreprocessing2(template,flag,mask,spp_json_directory,Result_directory,subFolderList,side,Sphingidae)
vdlist={'dorsal','ventral'};
maskf=mask;
%Find the symmetric axis
[symCentroid,symAxis,~]=findSymetricAxes(maskf);
disp('The symatric axis has been found.');
 %Find the center based on regionprop function   
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
 
 %Pick the most suitable one
 if cenDiff<50
    realCen= regioncen;
     disp('The centroid of the {entire mask} is used as the real centroid.');
 else
     realCen= cenregioncen;
      disp('The centroid of the {deleted central region} is used as the real centroid.');
 end
 disp('The centroid has been determined.');
 
%Derive the coordinates of corners
 ulCorner=boundingBox(1:2);
 urCorner=[boundingBox(1)+boundingBox(3),boundingBox(2)];
 llCorner=[boundingBox(1),boundingBox(2)+boundingBox(4)];
 lrCorner=[boundingBox(1)+boundingBox(3),boundingBox(2)+boundingBox(4)];
 allFrameCorners=[ ulCorner; urCorner; llCorner; lrCorner];
%%
%Prepare the symmetric axes for plotting
%create symmetric axes based on eigenvector
symOrtho=reshape(null(symAxis(:).'),1,[]);
dim_1=realCen+symAxis*size(maskf,1)/3;
dim_1plot=[realCen(1),dim_1(1);realCen(2),dim_1(2)];
dim_2=realCen+symOrtho*size(maskf,1)/9;
dim_2plot=[realCen(1),dim_2(1);realCen(2),dim_2(2)];

inspoutname=fullfile(Result_directory,'Shape_analysis',subFolderList{1},[template,'_',vdlist{side},flag,'_primary_key_pts.png']);
figinsp=figure('visible', 'off');
imshowpair(maskf,cenregion);
hold on;
plot(symCentroid(1),symCentroid(2),'y*');
plot(symCentroid(1),symCentroid(2),'bo');
%plot(symSkeletonCentroid(1),symSkeletonCentroid(2),'yo'); %Turn this on if the user wants to see this reference point
plot(regioncen(1),regioncen(2),'b*');
plot(cenregioncen(1), cenregioncen(2),'r+');
line(dim_1plot(1,:),dim_1plot(2,:), 'Color', 'r' ,'LineWidth', 1);
line(dim_2plot(1,:),dim_2plot(2,:), 'Color', 'b' ,'LineStyle','--','LineWidth', 1);
plot(allFrameCorners(:,1),allFrameCorners(:,2),'y+');
%rectangle('Position', boundingBox, 'EdgeColor','w','LineStyle','-.', 'LineWidth', 1); %Turn this on if the user wants to see this bounding box
%plot(realCen(1),realCen(2),'ro','LineWidth', 2); %Turn this on if the user wants to see the centroid
hold off;
export_fig(figinsp,inspoutname, '-png','-r150');
close(figinsp);
disp('An image indicating the primary key points of specimen image has been saved.');
%%
%Show the mask in an image
maskoutname=fullfile(Result_directory,'Shape_analysis',subFolderList{2},[template,'_',vdlist{side},flag,'_mask_cen.png']);
figmask=figure('visible', 'off');
imshow(maskf);
hold on;
plot(realCen(:,1),realCen(:,2),'bo','LineWidth', 1);
line(dim_1plot(1,:),dim_1plot(2,:), 'Color', 'b' ,'LineWidth', 1);
line(dim_2plot(1,:),dim_2plot(2,:), 'Color', 'b' ,'LineStyle','--','LineWidth', 1);
hold off;
export_fig(figmask,maskoutname, '-png','-r150');
close(figmask);
disp('An image indicating the mask and centroid of specimen image has been saved.');
%%
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
disp('########## Begin to find the corner between fore and hindwings. #########');
disp('Begin to find the corner between left fore and hindwings.');
nStrongCornersList=[500,1000,2000,3000,4000]; %shorter trail list can be applicable: [500,1000,2000,4000];
nSectionList=[60:-5:20]; %number of elements in the list should greater than 4  ([20:5:50];)

slopeSwitch='wingEdge';
[conjPt,forehindCorner,~]=findForeHindCorner(nStrongCornersList,nSectionList,maskf,realCen,symAxis,ulCornerDV,boundingBoxDV,slopeSwitch);
if length(forehindCorner(forehindCorner(:,1)>0))<5
    slopeSwitch='cenAxis';
    [conjPt,forehindCorner,~]=findForeHindCorner(nStrongCornersList,nSectionList,maskf,realCen,symAxis,ulCornerDV,boundingBoxDV,slopeSwitch);
end
forehindCornerL=conjPt;
disp('The corner between left fore and hindwings has been found.');
disp('Begin to find the corner between right fore and hindwings.');

slopeSwitch='wingEdge';
[conjPt,forehindCorner,~]=findForeHindCorner(nStrongCornersList,nSectionList,maskf,realCen,symAxis,lrCornerDV,boundingBoxDV,slopeSwitch);
if length(forehindCorner(forehindCorner(:,1)>0))<5
    slopeSwitch='cenAxis';
    [conjPt,forehindCorner,~]=findForeHindCorner(nStrongCornersList,nSectionList,maskf,realCen,symAxis,lrCornerDV,boundingBoxDV,slopeSwitch);
end
forehindCornerR=conjPt;
disp('The corner between right fore and hindwings has been found.');
disp('########## Two corners between fore and hindwings have been found. #########');
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
export_fig(figmask1,maskoutname1, '-png','-r150');
close(figmask1);
disp('An image indicating the mask and fore-hindwing corners has been saved.');
%%
%Find best beltWpar
disp('Start to find the optimal body width parameter.');
beltWpar=findbeltWpar(maskf,forehindCornerL,forehindCornerR,realCen,boundingBox);
disp('The body width parameter has been found.');
%%
disp('Start to crop upper and lower masks based on the fore-hindwing corners.');
upperSegMaskPts=[[0,0];[size(maskf,2),0];forehindCornerR;realCen;forehindCornerL];
upperSegMask = poly2mask(round(upperSegMaskPts(:,1)),round(upperSegMaskPts(:,2)),size(maskf,1),size(maskf,2));
upperMask=immultiply(maskf,upperSegMask);

lowerSegMaskPts=[[0,size(maskf,1)];[size(maskf,2),size(maskf,1)];forehindCornerR;realCen;forehindCornerL];
lowerSegMask = poly2mask(round(lowerSegMaskPts(:,1)),round(lowerSegMaskPts(:,2)),size(maskf,1),size(maskf,2));
lowerMask=immultiply(maskf,lowerSegMask);
disp('Upper and lower masks have been cropped out based on the fore-hindwing corners.');

%%
disp('########## Start to find the corner between wing and body. ##########');
disp('Begin to find the corner between the left forewing and the body.');
%Left forewing
nStrongCornersList=[1000,1500,2000,3000,4000]; %shorter trail list can be applicable: [1000,1500,2000,4000];
nSectionList=reshape( [[20:4:50];[50:-4:20]],1,[]); %use a mix parameter lists low-high + high-low to search for better results
[conjPt, conjCorners]=findBodyWingCorner2(nStrongCornersList,nSectionList,upperMask,realCen,symAxis,forehindCornerL,'LF',boundingBox,beltWpar,maskf);
conjCornerLF=conjPt;
disp('The corner between the left forewing and body has been found.');

disp('Begin to find the corner between the right forewing and the body.');
%Right forewing
[conjPt, conjCorners]=findBodyWingCorner2(nStrongCornersList,nSectionList,upperMask,realCen,symAxis,forehindCornerR,'RF',boundingBox,beltWpar,maskf);
conjCornerRF=conjPt;
disp('The corner between the right forewing and the body has been found.');

disp('Begin to find the corner between the left hindwing and the body.');
%Left hindwing
nSectionList=[20:4:40]; %number of elements in the list should greater than 4
[conjPt, conjCorners]=findBodyWingCorner2(nStrongCornersList,nSectionList,lowerMask,realCen,symAxis,forehindCornerL,'LH',boundingBox,beltWpar,maskf);
conjCornerLH=conjPt;
disp('The corner between the left hindwing and the body has been found.');

disp('Begin to find the corner between the right hindwing and the body.');
%Right hindwing
[conjPt, conjCorners]=findBodyWingCorner2(nStrongCornersList,nSectionList,lowerMask,realCen,symAxis,forehindCornerR,'RH',boundingBox,beltWpar,maskf);
conjCornerRH=conjPt;
disp('The corner between the right hindwing and the body has been found.');
disp('########## All landmark corners are found. ##########');

%Integrate the landmark points: (Left-ForeHind Corner, Left-Fore-Body corner, Right-Fore-Body corner, Right-Fore-Hind Corner,Right-Hind-Body corner, Left-Hind-Body corner)
segPts=[forehindCornerL;conjCornerLF;conjCornerRF;forehindCornerR;conjCornerRH;conjCornerLH];


%%
%Find the slopes of the wings' main axes: (LeftFore, LeftHind, RightFore, RightHind)
disp('########## Start to find the slope of long axes of wings. ##########');
[WingAxesSlopes,tipPts]=findSlopesOfWingAxes(maskf,symAxis,realCen,segPts,boundingBox);
WingAxesSlopes(isinf(WingAxesSlopes))=10^10; %a small trick to prevent Inf
disp('########## The slope of long axes of wings has been determined. ##########');
%%%%%%%%%%%%%%%Special Region For Sphingidae%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For Sphingidae use this approach to derive the wingAxesSlopes (no hindwing rotation)
%Sphingidae=1; %1: Sphingidae; 0: not Sphingidae
%This argument won't influence the final result. It only changes the
%visualization of preliminary analysis
if Sphingidae==1
    WingAxesSlopes(2)=0;
    WingAxesSlopes(4)=0;
    disp('The rotation of hindwings will not be corrected for Sphingidae.');
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
export_fig(figchar,charoutname, '-png','-r150');
close(figchar);
disp('An image indicating all key characters of the specimen image has been saved.');
%%
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
         disp('Complete with both masks and json file for segmentation');
    catch
        disp('Something went wrong in the analysis with json file, return to using naive segmentation');
        [wingParts,refPts]=segWings(mask,realCen,segPts);
    end
else
    disp('Use only mask for segmentation');
    [wingParts,refPts]=segWings(mask,realCen,segPts);
end
disp('########## Four wing parts are cropped. ##########');
%The following scripts are kept for visulization when debugging
%figure,imshow(mask); hold on; plot(segPts(:,1),segPts(:,2),'r+');
%figure,imshow(mask); hold on; plot(refPts(:,1),refPts(:,2),'r+');
disp('########## Start to rotate each wing part. ##########');
[tformWingImg, tformPtData]=wingRotation(wingParts,WingAxesSlopes,refPts);
disp('########## The rotation of all wings has been finished. ##########');

disp('########## Start to compose forewings and hindwings. ##########');
%rotate and normalize the position of a pair of forewings and hindwings
shapeImg=combineForeHindWings(tformWingImg,tformPtData);
disp('########## Forewings and hindwings have been composed. ##########');

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
export_fig(figshape,shapevisoutname, '-png','-r150');
close(figshape);
disp('An image showing adjusted specimen wing shape has been saved.');

shapevisoutname=fullfile(Result_directory,'Shape_analysis',subFolderList{4},[template,'_',vdlist{side},flag,'_wings_shape_n_pts.png']);
figshape=figure('visible', 'off');
imshow(NewAllComImg);
hold on;
plot(allinfo(:,1)-NewOrigion(1),allinfo(:,2)-NewOrigion(2),'yx','LineWidth', 2);
hold off;
export_fig(figshape,shapevisoutname, '-png','-r150');
close(figshape);
disp('An image showing adjusted specimen wing shape and boundary pts has been saved.');
end