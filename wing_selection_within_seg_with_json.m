function target=wing_selection_within_seg_with_json(mask,cropMask, trimPx, pick_criteria)
%     pick_criteria='upper';
%     pick_criteria='lower';

    
 if strcmp(pick_criteria,'upper')
    wingCandidates=bwareafilt(bwareaopen(logical(imdilate(imerode(immultiply(mask,cropMask),strel('disk',trimPx)),strel('disk',trimPx))),1000),2); %also eliminate small object smaller than 1000
elseif strcmp(pick_criteria,'lower')
    wingCandidates=bwareafilt(bwareaopen(logical(immultiply(mask,cropMask)),1000),2); %Preserve tail, so no erosion
else
    wingCandidates=bwareafilt(bwareaopen(logical(imdilate(imerode(immultiply(mask,cropMask),strel('disk',trimPx)),strel('disk',trimPx))),1000),2); %also eliminate small object smaller than 1000
end
    
    [bL,bn] = bwlabel(wingCandidates);
    stat=regionprops(bL, 'Centroid','Area');

    if bn==1
        target=wingCandidates;
    else
        centroid1=stat(1).Centroid;
        centroid2=stat(2).Centroid;
%         area1=stat(1).Area;
%         area2=stat(2).Area;

        if centroid1(2)>centroid2(2)
            upperID=2;
            lowerID=1;
        else
            upperID=1;
            lowerID=2;
        end

        if abs(centroid1(2)-centroid2(2))>50
            if strcmp(pick_criteria,'upper')
                target= bL==upperID;
            elseif strcmp(pick_criteria,'lower')
                target= bL==lowerID;
            else
                target= bwareafilt(wingCandidates,1);
            end
        else
            target= bwareafilt(wingCandidates,1);
        end
    end
end
