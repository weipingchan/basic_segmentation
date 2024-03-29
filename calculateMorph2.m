function calculateMorph2(spp_mat_directory,Code_directory,spp_json_directory,Result_directory,infilename,SphingidaeOrNot)
%The variable SphingidaeOrNot only influences the preliminary visualization,
%so it won't change the final result even when the user specifies 0 to a
%Sphingidae

if size(spp_mat_directory,2)==1 spp_mat_directory=spp_mat_directory{1};, end;
if size(Code_directory,2)==1 Code_directory=Code_directory{1};, end;
if size(Result_directory,2)==1 Result_directory=Result_directory{1};, end;
if size(infilename,2)==1 infilename=infilename{1};, end;
if ~isnumeric(SphingidaeOrNot) SphingidaeOrNot=str2num(SphingidaeOrNot);, end;

% Turn off this warning "Warning: Image is too big to fit on screen; displaying at 33% "
% To set the warning state, you must first know the message identifier for the one warning you want to enable. 
warning('off', 'Images:initSize:adjustingMag');

addpath(genpath(Code_directory)) %Add the library to the path

disp('Start to create / find corresponding folders.');
%Create result directory
if ~exist(fullfile(Result_directory,'Shape_analysis'), 'dir')
    mkdir(fullfile(Result_directory,'Shape_analysis'));
end

subFolderList={'primary_key_pts','mask','main_characters','wing_shape_visualization','shape_context','shape_mask','wing_shape_matrices-seg','wing_segmentation'};

for fold=1:length(subFolderList)
    if ~exist(fullfile(Result_directory,'Shape_analysis',subFolderList{fold}), 'dir')
        mkdir(fullfile(Result_directory,'Shape_analysis',subFolderList{fold}));
    end
    disp(['corresponding folder ', subFolderList{fold}, ' is created / found.']);
end

vdlist={'dorsal','ventral'};
disp('Start to read images into memory');
matinname=infilename;
if contains(matinname, 'dorsal')
    side=1;
    template0=strsplit(matinname,['_',vdlist{side}]);
    template=template0{1};

elseif contains(matinname, 'ventral')
    side=2;
    template0=strsplit(matinname,['_',vdlist{side}]);
    template=template0{1};
end

if contains(matinname, '_r_')
    flag='_r';
elseif contains(matinname, '_m_')
    flag='_m';
else
    flag='';
end

disp(['Start to analyze specimen: ', template,'_',vdlist{side},flag]);
sppmat0=load(fullfile(spp_mat_directory,matinname));
fieldName=cell2mat(fieldnames(sppmat0));
sppmat=sppmat0.(fieldName);
clear sppmat0;

mask=sppmat{end-1};
scale=sppmat{end};
clear sppmat;

[realCen, symAxis, symOrtho, boundingBox,  tipPts, refPts, wingParts, WingAxesSlopes, tformWingImg, shapeImg, allComImg, allinfo]=wingShapePreprocessing2(template,flag,mask,spp_json_directory,Result_directory,subFolderList,side,SphingidaeOrNot);

bodyTrimPx=5;
antTrimPx=10; %Adjust this to extract bold antennae (e.g. some moths') if necessary
[bodyMask, bodyCharacters, antennaMask,antennaCharacters]=body_antenna_module2(mask, wingParts, refPts, tipPts, bodyTrimPx, antTrimPx, scale);

%Generate a segmentation image for visualization
segmented_img=zeros(size(bodyMask));
segmented_img(mask==1)=7;
segmented_img(antennaMask==1)=6;
segmented_img(bodyMask==1)=5;
if side==1 %On dorsal side, forewings are dominant
    segmented_img(wingParts{2}==1)=3;
    segmented_img(wingParts{4}==1)=4;
    segmented_img(wingParts{1}==1)=1;
    segmented_img(wingParts{3}==1)=2;
elseif side==2 %On ventral side, hindwings are dominant
    segmented_img(wingParts{1}==1)=2;
    segmented_img(wingParts{3}==1)=1;
    segmented_img(wingParts{2}==1)=4;
    segmented_img(wingParts{4}==1)=3;
end 
%0 is background
%1 is left forewing from dorsal side
%2 is right forewing from dorsal side
%3 is left hindwing from dorsal side
%4 is right hindwing from dorsal side
%5 is body
%6 is antennae
%7 is uncertain

%Create visualization; reference: https://color.adobe.com/zh/explore
color_idx=[[0,[0,0,0]],
    [1,[242, 200, 5]],
    [2,[242, 135, 5]],
    [3,[242, 75, 153]],
    [4,[166, 3, 33]],
    [5,[242, 240, 240]],
    [6,[171, 5, 242]],
    [7,[40, 40, 40]]
    ];

segmented_img_vis=label_segmentation_visualization(segmented_img,color_idx);

%Generate output matrix
morphinfo=cell(0,15);
morphinfo{1}=mask; %Original mask
morphinfo{2}=realCen; %Original centroid.
morphinfo{3}=symAxis; %The vertical symmetric axis
morphinfo{4}=symOrtho; %The horizontal symmetric axis
morphinfo{5}=tipPts; %tips of forewings
morphinfo{6}=refPts; %Important segment points
morphinfo{7}=wingParts; %right and left X fore and hindwings. Total: 4
morphinfo{8}=WingAxesSlopes; %Slopes for the front edge of forewings and rear edge of hindwings
morphinfo{9}=tformWingImg; %Re-posed wings
morphinfo{10}=shapeImg; %Re-posed fore-hindwings in the same panel and related points in new coordinations.
morphinfo{11}=allinfo; %Reference points before and after realignment
morphinfo{12}=scale; %The scale bar (number of pixels = 1 cm)
morphinfo{13}=segmented_img; %The segmentation image for visualization
morphinfo{14}=bodyCharacters; %The body length and width (in cm)
morphinfo{15}=antennaCharacters; %Antennae length, width, bolb width, degree of curvature (all in millimeters); first row is left one, second is right one.

matoutname=fullfile(Result_directory,'Shape_analysis',subFolderList{7},[template,'_',vdlist{side},flag,'_morph-seg.mat']);
save(matoutname,'morphinfo'); %save the specimen matrix

%Save the segmented image
segvisoutname=fullfile(Result_directory,'Shape_analysis',subFolderList{8},[template,'_',vdlist{side},flag,'_wing_segmentation.jpg']);
figseg=figure('visible', 'off');
imshow(segmented_img_vis)
export_fig(figseg,segvisoutname, '-jpg','-r150');
close(figseg);
disp('An image showing img segmentation has been saved.');

%Move those images that have been analyzed to a subdirectory
finishedDir='done';
if ~exist(fullfile(spp_mat_directory,finishedDir), 'dir')
    mkdir(fullfile(spp_mat_directory,finishedDir));
end
movefile(fullfile(spp_mat_directory,matinname),fullfile(spp_mat_directory,finishedDir));
disp('Analyzed images have been moved to done directory.');

end