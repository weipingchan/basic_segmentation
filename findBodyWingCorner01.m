function [conjPt, conjCorners]=findBodyWingCorner01(nStrongCornersList,nSectionList,partMask,realCen,symAxis,tarCorner,LeftRightForeHind,boundingBox,beltWpar,maskf)
disp('##############################'); 
disp('Run [find_body_wing_corner3] module');
    nsecn=1;
    conjCorners=zeros(length(nSectionList),2);
    for nstrongID=1:length(nStrongCornersList)
        nStrongCorners=nStrongCornersList(nstrongID);
        disp(['Try number of interesting points: ',num2str(nStrongCorners)]);
        for nsec=1: length(nSectionList)
            nSection=nSectionList(nsec);
            try
                disp(['Try section number: ',num2str(nSection)]);
                conjCorners0=find_body_wing_corner3(partMask,nStrongCorners,realCen,symAxis,tarCorner,LeftRightForeHind,nSection,boundingBox,beltWpar);
                if ~isempty(conjCorners0)
                    conjCorners(nsecn,:)=conjCorners0;
                    nsecn=nsecn+1;
                end
            catch
                disp('DOESNOT work.');
                %conjCorners(nsecn,:)=[0 0];
                %break
            end
        end
        if nnz(conjCorners(:,1))>=13
            disp([num2str(nStrongCorners),' interesting points WORKS.']);
            disp(['Parameter [nStrongCorners]: ',num2str(nStrongCorners)]);
            break
        else
            disp([num2str(nStrongCorners),' interesting points does not work.']);
        end
    end
    
    %%
    %if there are two set of points having the same frequency, use the
    %one near the tarCorner
    if LeftRightForeHind=='LF'
        comCorner=[1, 1];
    elseif LeftRightForeHind=='RF'
        comCorner=[size(partMask,2), 1];
    elseif LeftRightForeHind=='LH'
        comCorner=[1,size(partMask,1)];
    else
        comCorner=size(partMask);
    end
    %%
    %Search among top 2 highly-frequent candidate points to determine the
    %final one
    [Au,~,ic] = unique(conjCorners(conjCorners(:,1)>0,:),'rows');
    countAu = accumarray(ic,1);
    
    if isempty(Au)
        %%Provide a temperary point
        %Calculate necessary vectors
        partMask2=bwareafilt(logical(imdilate(imerode(partMask,strel('disk',5)),strel('disk',5))),[100,Inf]); %remove antanee
        verVector=symAxis*size(partMask2,1);
          
        %Find all edge points based on de antanee mask
        [specimenB,~]=bwboundaries(partMask2);

        sppEdgePt=[];
        for tid=1:length(specimenB)
            sppEdgePt=[sppEdgePt; specimenB{tid}];
        end

        if contains(LeftRightForeHind, 'F')
            newCen=realCen+(tarCorner-realCen)/4+[0,-5]; %find a new starting centroid point
            tmpSegPtsf=[newCen-verVector ; newCen+verVector];
            [intersectXf,intersectYf]= polyxpoly(tmpSegPtsf(:,1),tmpSegPtsf(:,2),sppEdgePt(:,2),sppEdgePt(:,1));
            seg2f=[intersectXf,intersectYf];
            seg2f=seg2f(seg2f(:,2)<newCen(2),:);
        else
            newCen=realCen+(tarCorner-realCen)/4+[0,5]; %find a new starting centroid point
            tmpSegPtsf=[newCen-verVector ; newCen+verVector];
            [intersectXf,intersectYf]= polyxpoly(tmpSegPtsf(:,1),tmpSegPtsf(:,2),sppEdgePt(:,2),sppEdgePt(:,1));
            seg2f=[intersectXf,intersectYf];
            seg2f=seg2f(seg2f(:,2)>newCen(2),:);
        end
        travelDist = sqrt(sum(bsxfun(@minus,seg2f,newCen).^2,2));
    
        conjPt=seg2f(travelDist==min(travelDist),:);
        disp('Assign a key point.');
%         figure,imshow(partMask);hold on;
%         plot(conjPt(:,1),conjPt(:,2),'rO');
%         plot(newCen(:,1),newCen(:,2),'bO');
%         plot(realCen(:,1),realCen(:,2),'b*');
%         plot(tmpSegPtsf(:,1),tmpSegPtsf(:,2),'r');
% %         plot(sppEdgePt(:,2),sppEdgePt(:,1),'y*')
    else    
        %Find the top 2 value without 
        countAu0=sort(countAu,'descend');
        if length(countAu0)>=2
            commonest=Au(countAu>1&(countAu==countAu0(1)|countAu==countAu0(2)),:); %Chose top 2 highly-frequent pts (frequency >1)
        else
            try
                commonest=Au(countAu>1&countAu==countAu0(1),:); %if there is only one highly-frequent pt (frequency >1)
            catch
                commonest=[]; %This line is added to prevent interuption of running
            end
        end

        if size(commonest,1)>1
            disp(['Find ',num2str(size(commonest,1)),' candidate points; begin refining.']);
            if (countAu0(1)-countAu0(2))/countAu0(1)>=0.5
                conjPt=Au(countAu==countAu0(1),:);
            else
                conjPt=findCloestPt(commonest,comCorner);
            end
        elseif isempty(commonest)
            conjPt=Au(round(size(Au,1)/2),:);
        else
            conjPt=commonest;
        end
        disp('Find the key point.');
    end
end