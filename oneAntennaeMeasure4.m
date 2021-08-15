function [bolbMorph,antFig]=oneAntennaeMeasure4(mask, forkPt00, oneAntenae,antenna1,It0,Dd0,realCen,scalelen)
     try %prevent protential error and interoption
        It = bwmorph(antenna1,'thin','inf');
        [Bbx,Bby] =  find(bwmorph(It,'branchpoints'));
        %select branch point and convert to index
        dist = sqrt(sum(bsxfun(@minus, [Bbx,Bby], flip(forkPt00)).^2,2));
        [~,minLoc]=min(abs(dist-mean(dist(dist<mean(dist)+std(dist)*1/2,:)))); %find the one closest to the mean of those points near forkPt
        Bbi=sub2ind(size(antenna1),Bbx(minLoc),Bby(minLoc)); %The forkPoint we want

        [i,j] = find(immultiply(bwmorph(It,'endpoints'),oneAntenae)); %Calculation with refined end points
        Dd = bwdistgeodesic(It,Bbi,'quasi');
        distanceAnt=[];
        for n = 1:numel(i)
            distanceAnt=[distanceAnt; [i(n),j(n),Dd(i(n),j(n))]];
        end
        %%
        if size(distanceAnt,1)>1
            %Cluster points having long distance to the fork
            Dist = sqrt(sum(bsxfun(@minus, distanceAnt(:,1:2), flip(forkPt00)).^2,2));
            antennaTipDist0= sortrows(distanceAnt(Dist>mean(Dist)+std(Dist)*1/2,:),2); %first one is the left one
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
        else
            antennaTipDist=distanceAnt;
        end
    catch
        antennaTipDist=[];
    end


    if ~isempty(antennaTipDist)
        try
            antenaeEdgeMask=imdilate(oneAntenae,strel('disk',1))-imerode(oneAntenae,strel('disk',1));
            [ei,ej] = find(immultiply(It0,antenaeEdgeMask));
            antennaeBase0=sortrows([ei,ej],-1); %sort in descending order
            antennaeBase=antennaeBase0(1,:);
            fork2baseDist=Dd0(antennaeBase(1),antennaeBase(2));
        catch
            fork2baseDist=0;
            antennaeBase=forkPt00;
        end
        antL= antennaTipDist(1,3)-fork2baseDist;
    else
        antL= [];
        antennaeBase=[];
    end
    
    try
        [bolbMorph1,~,tipAntenae,tipAnt]=oneAntennaeMeasure03(mask, forkPt00, oneAntenae,antL, antennaeBase);
    catch
        bolbMorph1=-9999 + zeros(1, 4);
        tipAnt=[];
    end
    %judge if the antennae is the right side or left side
    ss01 = regionprops(oneAntenae,'centroid');
    antCens=ss01.Centroid;
    if antCens(1)<realCen(1)
        rfAnt=1;
    else
        rfAnt=2;
    end
    bolbMorph=-9999 + zeros(2,4);
    if ~isempty(antL) && ~isempty(tipAnt)
        basetipdist=pdist([tipAnt;antennaeBase],'euclidean');
        curveDegree=antL/basetipdist;
    elseif isempty(antL) && ~isempty(tipAnt) %Added March 14, 2020
        if bolbMorph1(1,1)<0
            antL=-9999;
            curveDegree=-9999;
        else
            antL=bolbMorph1(1,1);
            curveDegree=bolbMorph1(1,4);
        end
    else  %Added March 14, 2020
        antL=-9999; %Added March 14, 2020
        curveDegree=-9999; %Added March 14, 2020
    end
    if antL>0 antLout=antL/scalelen*10;, else antLout=-9999;,  end;
    if bolbMorph1(2)>0 antWout=bolbMorph1(2:3)/scalelen*10;, else antWout=bolbMorph1(2:3);,  end;
    bolbMorph(rfAnt,:)=[antLout, antWout, curveDegree];
    antFig=oneAntenae;
            
            %%
%         %Plot antenna
%         box=regionprops(oneAntenae,'BoundingBox'); %for visualization
%         antBasePts=antennaeBase; %for visualization
%         textPos=box.BoundingBox(1:2); %for visualization
%         
%         antFigv=figure;
%         imshow(labeloverlay(double(mask*0.1+oneAntenae*0.5),It,'Colormap','autumn','Transparency',0));hold on;
% %     plot(antBasePts(:,2), antBasePts(:,1),'yO');
%         plot(antBasePts(:,2), antBasePts(:,1),'y*');
%             if rfAnt==1
%                 text(double(textPos(1))-30,double(textPos(2))-10,[num2str(round(bolbMorph(rfAnt,1),2)),'mm'],'color','g');
%             else
%                 text(double(textPos(1))+30,double(textPos(2))-10,[num2str(round(bolbMorph(rfAnt,1),2)),'mm'],'color','g');
%             end
%         hold off;   
%         
%         %%
%         %Plot tips
%         box=regionprops(oneAntenae,'BoundingBox'); %for visualization
%         antBasePts=antennaeBase; %for visualization
%         textPos=box.BoundingBox(1:2); %for visualization
%         
%         antFigv=figure;
%         imshow(labeloverlay(double(mask*0.1+oneAntenae*0.2+tipAntenae),It,'Colormap','autumn','Transparency',0));hold on;
% %     plot(antBasePts(:,2), antBasePts(:,1),'yO');
%         plot(antBasePts(:,2), antBasePts(:,1),'yx');
%             if rfAnt==1
%                 text(double(textPos(1))-30,double(textPos(2))-10,[num2str(round(bolbMorph(rfAnt,3),2)),'mm'],'color','w');
%             else
%                 text(double(textPos(1))+30,double(textPos(2))-10,[num2str(round(bolbMorph(rfAnt,3),2)),'mm'],'color','w');
%             end
%         hold off;   
end