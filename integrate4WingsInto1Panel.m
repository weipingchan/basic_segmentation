function [allComImg, allinfo]=integrate4WingsInto1Panel(shapeImg,tformWingImg,refPts)
%Fore wing mix
disp('Start to combine fore wings.');
Lpar=shapeImg{1,2};
Rpar=shapeImg{2,2};
bodyW=pdist([refPts(end-2,:);refPts(end,:)],'euclidean')+200; %This value determine the distance between left and right wings
Lref=Lpar(2,:)+[bodyW 0];
Rref=Rpar(2,:);
disp('The reference image is Left Fore wing.');

LforImg=tformWingImg{1};
RforImg=tformWingImg{3};

horMvDis=Lref-Rref;
if horMvDis(1)<0 || horMvDis(2)<0 %If the moving direction is negetive, keping the full size will chnge the coordination of the target object
    mvRImg = imtranslate(RforImg,horMvDis);
    disp('Full output view is NOT applied since variable [horMvDis] containing negative value.');
else
    mvRImg = imtranslate(RforImg,horMvDis,'FillValues',0,'OutputView','full');
end
disp('Right Fore wing has been re-registered.');
FcomImg = imfuse(LforImg,mvRImg,'blend');
disp('The composition of fore wings is finished.');

%Hind wing mix
disp('Start to combine hind wings.');
LhindImg=tformWingImg{2};
RhindImg=tformWingImg{4};

LverMvDis=Lpar(2,:)-(Lpar(3,:));
if LverMvDis(1)<0 || LverMvDis(2)<0
    LmvRImg = imtranslate(LhindImg,LverMvDis);
    disp('Full output view is NOT applied since variable [LverMvDis] containing negative value.');
else
    LmvRImg = imtranslate(LhindImg,LverMvDis,'OutputView','full');
end
disp('Left Hind wing has been re-registered.');

RverMvDis=Rpar(2,:)-Rpar(3,:);
if RverMvDis(1)<0 || RverMvDis(2)<0
    RmvRImg0 = imtranslate(RhindImg,RverMvDis);
    disp('Full output view is NOT applied since variable [RverMvDis] containing negative value.');
else
    RmvRImg0 = imtranslate(RhindImg,RverMvDis,'OutputView','full');
end
if horMvDis(1)<0 || horMvDis(2)<0
    RmvRImg = imtranslate(RmvRImg0,horMvDis);
    disp('Full output view is NOT applied since variable [horMvDis] containing negative value.');
else
    RmvRImg = imtranslate(RmvRImg0,horMvDis,'OutputView','full');
end
disp('Right Hind wing has been re-registered.');

HcomImg = imfuse(LmvRImg,RmvRImg,'blend');
disp('The composition of hind wings is finished.');
%Combine both image to derive the final one
allComImg = imfuse(FcomImg,HcomImg);
disp('The composition of all wings is finished.');

%derive new Info
allinfo=[Lpar(1,:) ; Lpar(4,:) ; Rpar(1,:)+horMvDis; Rpar(4,:)+horMvDis];
disp('The composition of all parameters is generated.');
end