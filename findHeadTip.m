function headTipPt=findHeadTip(body)
%find the head instead of fork point
        [sideE0,~]=bwboundaries(body);
        sideEdgePt=sideE0{1};
        
        [sideRecx,sideRecy,~,~] = minboundrect(sideEdgePt(:,1),sideEdgePt(:,2));

        recSlope=round((diff(sideRecy) ./ diff(sideRecx)),4);
        candidateLoc=find(abs(recSlope)==max(abs(recSlope)));
        
        Recxy= [sideRecx,sideRecy];
        horLa=Recxy(candidateLoc(1):candidateLoc(1)+1,:);
        horLb=Recxy(candidateLoc(2):candidateLoc(2)+1,:);
        
        if mean(horLa(:,1)) >mean(horLb(:,1))
            topL=horLb;
            botL=horLa;
        else
            topL=horLa;
            botL=horLb;
        end
                
        edgeDist = [];
        for p=1:length(sideEdgePt)
            pd = point_to_line(sideEdgePt(p,:), topL(1,:), topL(2,:));
            edgeDist = [edgeDist; pd];
        end
        
        %find the smallest distance and use that as an index into B:
        cloestRecPt0=flip(sideEdgePt(edgeDist==min(edgeDist),:),2);
        if size(cloestRecPt0,1)>1
            headTipPt=cloestRecPt0(1,:);
        else
            headTipPt=cloestRecPt0;
        end

end