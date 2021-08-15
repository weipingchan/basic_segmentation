function [bolbMorph,antFig]=antennaMeasure8(mask,refPts,tipPts,wingMask,body, minimalAntennaLength, bodyraw, upperSegMask, scalelen)
    %The unit for all antenna parameters is mm
    % minimalAntennaLength=80; %The length of antenna should longer than this length or will be neglected

    realCen=refPts(8,:);
    forehindCornerL=tipPts(1,:);
    forehindCornerR=tipPts(2,:);

    % upperSegMaskPts=[[0,0];[size(mask,2),0];[size(mask,2),forehindCornerR(2)];forehindCornerR;realCen;forehindCornerL;[0,forehindCornerL(2)]];
    % upperSegMask = poly2mask(round(upperSegMaskPts(:,1)),round(upperSegMaskPts(:,2)),size(mask,1),size(mask,2));
    % antennaraw=mask-body-wingMask;
    % antenna=bwareafilt(logical(immultiply(antennaraw,upperSegMask)),1);

    antenna0=bwareafilt((mask-wingMask)>0,1);
    BWant0 = bwskel(logical(antenna0),'MinBranchLength',30);
    branchimage0 = bwmorph(BWant0, 'branchpoints');
    ss = regionprops(branchimage0,'centroid');
    forkPts0=cat(1,ss.Centroid);

    % figure,imshowpair(antenna0,BWant0);hold on;
    % plot(headTipPt(1),headTipPt(2),'rx');
    % plot(forkPt00(1),forkPt00(2),'ro');
    % plot(forkPts0(:,1),forkPts0(:,2),'b+');

    if ~isempty(forkPts0)
        %Find the closet point to the tip of the body
        try
            headTipPt=findHeadTip(body);
        catch
            headTipPt=[mean(forkPts0(:,1)),min(forkPts0(:,2))];
        end
        forkPts00=forkPts0(forkPts0(:,2)>headTipPt(:,2),:);
        lfLoc=forkPts00(:,1);
        midlineDist=abs(lfLoc-mean(lfLoc(forkPts00(:,2)>=median(forkPts00(:,2)))));
        forkPts01=forkPts00(midlineDist<=mean(midlineDist),:);
        forkPt00=findCloestPt(forkPts01,headTipPt);
    %     forkPt00=findCloestPt(forkPts0,headTipPt);
    %     [~, Locs] = ismember(forkPt00, forkPts0, 'rows');
    %     forkPt=findCloestPt(forkPts0,[mean(forkPts0(:,1)),min(forkPts0(:,2))]);   %verison 2
    %     forkPt=forkPts0(forkPts0(:,2)==min(forkPts0(:,2)),:);   %verison 1
    else
        forkPt00=flip(realCen);
    end

    extendingFactor=50;
    antennaMaskPts=[[0,0];[size(mask,2),0];[size(mask,2),forehindCornerR(2)];forehindCornerR;forkPt00+[0,extendingFactor];forehindCornerL;[0,forehindCornerL(2)]];
    antennaMask = poly2mask(round(antennaMaskPts(:,1)),round(antennaMaskPts(:,2)),size(mask,1),size(mask,2));

    antenna1=bwareafilt(logical(immultiply(antenna0,antennaMask)),1);

    %extract antenna image
    antenna02=bwareafilt(bwareaopen(immultiply(antenna1-body,antennaMask),200),[5, Inf],4); %Remove 1 or 2 pixels attachment
    cc0 = bwconncomp(antenna02,4);
    numberOfAntenna0  = cc0.NumObjects;
   
    try %prevent protential error and interoption
        It0 = bwmorph(antenna1,'thin','inf');
        [Bbx,Bby] =  find(bwmorph(It0,'branchpoints'));
        %select branch point and convert to index
        dist = sqrt(sum(bsxfun(@minus, [Bbx,Bby], flip(forkPt00)).^2,2));
        [~,minLoc]=min(abs(dist-mean(dist(dist<mean(dist)+std(dist)*1/2,:)))); %find the one closest to the mean of those points near forkPt
        Bbi=sub2ind(size(antenna1),Bbx(minLoc),Bby(minLoc)); %The forkPoint we want

        [i,j] = find(bwmorph(It0,'endpoints'));
        Dd0 = bwdistgeodesic(It0,Bbi,'quasi');
        distanceAnt0=[];
        for n = 1:numel(i)
            distanceAnt0=[distanceAnt0; [i(n),j(n),Dd0(i(n),j(n))]];
        end
        %%
        %Cluster points having long distance to the fork
        Dist = sqrt(sum(bsxfun(@minus, distanceAnt0(:,1:2), flip(forkPt00)).^2,2));
        antennaTipDist0= sortrows(distanceAnt0(Dist>mean(Dist)+std(Dist)*1/2,:),2); %first one is the left one
        antennaTipDistA0=sortrows(antennaTipDist0(antennaTipDist0(:,1)<mean(antennaTipDist0(:,1)),:),3);
        antennaTipDistB0=sortrows(antennaTipDist0(antennaTipDist0(:,1)>=mean(antennaTipDist0(:,1)),:),3);

        if ~isempty(antennaTipDistA0)
            antennaTipDistA=antennaTipDistA0(end,:);
        else
            antennaTipDistA=antennaTipDistA0;
        end
        if ~isempty(antennaTipDistB0)
            antennaTipDistB=antennaTipDistB0(end,:);
        else
            antennaTipDistB=antennaTipDistB0;
        end
        antennaTipDist=[antennaTipDistA;antennaTipDistB];
    catch
        antennaTipDist=[];
    end
    
    
%%
if ~isempty(antennaTipDist) && numberOfAntenna0~=0
    numberOfAntenna=size(antennaTipDist,1);
        
    if max(antennaTipDist(:,3))>minimalAntennaLength %Define antenna as those longer than minimalAntennaLength pixels
        if numberOfAntenna==2 && numberOfAntenna0>=2 %two antenae
            disp('TWO antenna are found');
            [bolbMorph,antFig]=twoAntennaMeasure5(mask, body, antenna1,distanceAnt0,forkPt00, bodyraw, upperSegMask, scalelen);
        else %Only one antenae
            disp('ONE anteanne is found');
            antennaraw=mask-body-wingMask;
            try
                oneAntenae=bwareafilt(logical(immultiply(antennaraw,antenna1)),1);
            catch
                oneAntenae=bwareafilt(logical(antennaraw),1);
            end
            try
                [bolbMorph,antFig]=oneAntennaeMeasure4(mask, forkPt00, oneAntenae,antenna1,It0, Dd0,realCen,scalelen);
           catch
                disp('Ops! There is actually NO antennae found.');
                bolbMorph=-9999 + zeros(2, 4);
                antFig=double(mask*0.1);
            end
        end

    else %No antenna
        bolbMorph=-9999 + zeros(2, 4);
        antFig=double(mask*0.1+antenna1*0.5);
    %     %Plot antenna
    %     antFig=figure('visible', 'off');
    %     imshow(labeloverlay(double(mask*0.1+antenna1*0.5),It,'Colormap','autumn','Transparency',0));
    %     hold off;
        disp('No antenna');
    end
else
    try
        antennaraw=mask-body-wingMask;
        oneAntenae=bwareafilt(logical(immultiply(antennaraw,antennaMask)),1);
        [bolbMorph,antFig]=oneAntennaeMeasure4(mask, forkPt00, oneAntenae,antenna1,It0, Dd0,realCen,scalelen);
        disp('ONE anteanne is found');
    catch
        bolbMorph=-9999 + zeros(2, 4);
        antFig=double(mask*0.1);
    %     %Plot antenna
    %     antFig=figure('visible', 'off');
    %     imshow(mask*0.1);
    %     hold off;
    disp('No antenna');
    end
end

end