function [bolbMorph,antFig,tipAntenae,tipAnt]=oneAntennaeMeasure03(mask, forkPt00, oneAntenae,antL0,antennaeBase0)
%A funciton to measure the morphology of one antenna

        antennae = oneAntenae;
        %%
        forkPt=forkPt00;
        
        %%
        It = bwmorph(antennae,'thin','inf');
        Bb =  bwmorph(It,'branchpoints');
               
        [i,j] = find(bwmorph(It,'endpoints'));
        Dd = bwdistgeodesic(It,find(Bb),'quasi');
        distanceAnt=[];
        for n = 1:numel(i)
            distanceAnt=[distanceAnt; [i(n),j(n),Dd(i(n),j(n))]];
        end
        Dist = sqrt(sum(bsxfun(@minus, distanceAnt(:,1:2), flip(forkPt)).^2,2));
        scalelist=[2/3,3/5,1/2,2/5];
        scaleID=1;
        while 1
            antennaTipDist2= sortrows(distanceAnt(Dist>mean(Dist)+std(Dist)*scalelist(scaleID),:),2);
            if ~isempty(antennaTipDist2) || scaleID>=length(scalelist)
                break
            else
                scaleID=scaleID+1;
            end
        end
        tipAnt= antennaTipDist2(1,1:2);
        if ~isempty(antL0) %Set up the condition flag
            if antL0>0
                antLflag=1;
            else
                antLflag=0;
            end
        else
            antLflag=0;
        end

        if ~isempty(antennaeBase0) %Set up the condition flag
            antBflag=1;
        else
            antBflag=0;
        end
        
        if antLflag==1 && antBflag==1
            antL= antL0;
            antennaeBase=antennaeBase0;
        else
            baseCandidate=distanceAnt(:,1:2);
            [Lia, Locb] = ismember(tipAnt, baseCandidate, 'rows');
            baseCandidate(Locb,:)=[]; %Remove the row
            Dd2 = bwdistgeodesic(It,tipAnt(2),tipAnt(1),'quasi');
            distanceAnt2=[];
            for n = 1:size(baseCandidate,1)
                distanceAnt2=[distanceAnt2; [baseCandidate(n,1),baseCandidate(n,2),Dd2(baseCandidate(n,1),baseCandidate(n,2))]];
            end
            antL= min(distanceAnt2(:,3));
            antennaeBase=distanceAnt2(distanceAnt2(:,3)==antL,1:2);
        end
        
        antB=regionprops(antennae,'Area','Perimeter');
        meanW= antB.Area/antL;
        basetipdist=pdist([tipAnt;antennaeBase],'euclidean');
        curveDegree=antL/basetipdist;        
        tipMask = createCirclesMask(mask, flip(tipAnt), antL/5);     %Create tip mask again (use 1/5 of total length), and use regionProp find width
        tipAntenae=bwareafilt(logical(immultiply(antennae,tipMask)),1);
        bolbB=regionprops(tipAntenae,'MinorAxisLength');
        bolbW=bolbB.MinorAxisLength;
        bolbMorph=[antL,meanW, bolbW, curveDegree];
        antFig=double(mask*0.1+antennae*0.8); 
        
%         %Plot antenna; script is preserved for visualization when debugging
%         antFig=figure('visible', 'off');
%         imshow(labeloverlay(double(mask*0.1+antennae*0.5),It,'Colormap','autumn','Transparency',0));hold on;
%         plot( forkPt2(:,1), forkPt(:,2),'bO');
%         plot( forkPt2(:,1), forkPt(:,2),'bx');
%         for n = 1:size(antennaTipDist2,1)
%             text(double(antennaTipDist2(n,2)),double(antennaTipDist2(n,1)),num2str(antennaTipDist2(n,3)),'color','g');
%         end
%         hold off;
end