function [tformWingImg, tformPtData]=wingRotation(wingParts,WingAxesSlopes,refPts)
wingPartNameList={'Left Fore Wing','Left Hind Wing','Right Fore Wing','Right Hind Wing'};
    tformWingImg=cell(size(wingParts,2),0);
    tformPtData=cell(size(wingParts,2),0);
    for pInd=1:4
        partSlope=WingAxesSlopes(pInd);
        partWing=wingParts{pInd};
        disp(['Start to rotate ',wingPartNameList{pInd},'. pInd: ',num2str(pInd)]);
        if pInd==1||pInd==3 %Fore wings
            if pInd==1 %Left fore wing
                partVector=[-1, -partSlope];
                hor = [-1, 0];
            elseif pInd==3 %Right fore wing
                partVector=[1, partSlope];
                hor = [1, 0];
            end
            %refVector=hor;
            [rpartWing, newRefPts]=rotateWings(pInd,partWing,partVector,hor,refPts);
            tformWingImg{pInd}=rpartWing;
            tformPtData{pInd}=newRefPts;
        else
            %Hind wings
            if partSlope~=0
                partVector=[sign(partSlope), abs(partSlope)];
                ver = [0, 1];
                [rpartWing, newRefPts]=rotateWings(pInd,partWing,partVector,ver,refPts);
            else %No rotation for hind wing (Eg. Sphingidae)
                disp('No rotation applied to hind wings');
                rpartWing=partWing;
                newRefPts=refPts;
            end
            %refVector=ver;
            tformWingImg{pInd}=rpartWing;
            tformPtData{pInd}=newRefPts;
        end
    end
end