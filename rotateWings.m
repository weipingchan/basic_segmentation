function [rpartWing, newRefPts]=rotateWings(pInd,partWing,partVector,refVector,refPts)
%A function to rotate wing pieces and the corresponding coordinates of key
%points; this is a function embedded in another function named wingRotation
    CosTheta = dot(partVector,refVector)/(norm(partVector)*norm(refVector));
    if pInd==1||pInd==3
        ThetaInDegrees = sign(partVector(1)*partVector(2))*acosd(CosTheta);
    else
        ThetaInDegrees = -sign(partVector(1)*partVector(2))*acosd(CosTheta);
    end
    tform = affine2d([cosd(ThetaInDegrees) -sind(ThetaInDegrees) 0; sind(ThetaInDegrees) cosd(ThetaInDegrees) 0; 0 0 1]);
    [rpartWing, partRef] = imwarp(partWing,tform);

    %transform the coordinates of key points as well
    [x1tr,y1tr]=transformPointsForward(tform,refPts(:,1),refPts(:,2));
    newRefPts=zeros(size(refPts,1),2);
    newRefPts(:,1) = x1tr - partRef.XWorldLimits(1);
    newRefPts(:,2) = y1tr - partRef.YWorldLimits(1);
end