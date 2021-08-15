function [rpartWing, newRefPts]=rotateWings(pInd,partWing,partVector,refVector,refPts)
    CosTheta = dot(partVector,refVector)/(norm(partVector)*norm(refVector));
    if pInd==1||pInd==3
        ThetaInDegrees = sign(partVector(1)*partVector(2))*acosd(CosTheta);
    else
        ThetaInDegrees = -sign(partVector(1)*partVector(2))*acosd(CosTheta);
    end
    %disp(['Variable [CosTheta]: ',num2str(CosTheta)]);
    %disp(['Variable [ThetaInDegrees]: ',num2str(ThetaInDegrees)]);
    tform = affine2d([cosd(ThetaInDegrees) -sind(ThetaInDegrees) 0; sind(ThetaInDegrees) cosd(ThetaInDegrees) 0; 0 0 1]);
    [rpartWing, partRef] = imwarp(partWing,tform);

    %transform the important dataset too
    [x1tr,y1tr]=transformPointsForward(tform,refPts(:,1),refPts(:,2));
    newRefPts=zeros(size(refPts,1),2);
    newRefPts(:,1) = x1tr - partRef.XWorldLimits(1);
    newRefPts(:,2) = y1tr - partRef.YWorldLimits(1);
    %disp(['The corresponding rotation of important landmark points have been done. ']);
end