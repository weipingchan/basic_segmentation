function [tformWingImg, tformPtData]=wingRotation(wingParts,WingAxesSlopes,refPts)
%A function to rotate each piece of wing and the corresponding coordinates of key points
wingPartNameList={'Left Fore Wing','Left Hind Wing','Right Fore Wing','Right Hind Wing'};
    tformWingImg=cell(size(wingParts,2),0);
    tformPtData=cell(size(wingParts,2),0);
    for pInd=1:4
        partSlope=WingAxesSlopes(pInd);
        partWing=wingParts{pInd};
        disp(['Start to rotate ',wingPartNameList{pInd},'. pInd: ',num2str(pInd)]);
        if pInd==1||pInd==3 %Forewings
            if pInd==1 %Left forewing
                partVector=[-1, -partSlope];
                hor = [-1, 0];
            elseif pInd==3 %Right forewing
                partVector=[1, partSlope];
                hor = [1, 0];
            end
            %refVector=hor;
            [rpartWing, newRefPts]=rotateWings(pInd,partWing,partVector,hor,refPts);
            tformWingImg{pInd}=rpartWing;
            tformPtData{pInd}=newRefPts;
        else
            %Hindwings
            if partSlope~=0
                partVector=[sign(partSlope), abs(partSlope)];
                ver = [0, 1];
                [rpartWing, newRefPts]=rotateWings(pInd,partWing,partVector,ver,refPts);
            else %No rotation for hindwing (Eg. Sphingidae)
                disp('No rotation applied to hindwings');
                rpartWing=partWing;
                newRefPts=refPts;
            end
            %refVector=ver;
            tformWingImg{pInd}=rpartWing;
            tformPtData{pInd}=newRefPts;
        end
    end
end