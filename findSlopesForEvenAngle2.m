function Ls2cen=findSlopesForEvenAngle2(symAxis,uSlopeL,nSection,LeftRightForeHind)

Ssym=-symAxis(2)/symAxis(1);
Sslo=-uSlopeL;
if LeftRightForeHind=='LF'
    symAxisTheta=mod(atan(Ssym),pi);
    uSlopeLTheta=atan(Sslo);
    Ls2cenTheta= linspace(symAxisTheta,uSlopeLTheta,nSection);
elseif LeftRightForeHind=='RF'
    symAxisTheta=mod(atan(Ssym),pi);
    if Sslo<0
        uSlopeLTheta=mod(atan(Sslo),pi);
    else
        uSlopeLTheta=atan(Sslo)+pi;
    end
    Ls2cenTheta= linspace(symAxisTheta,uSlopeLTheta,nSection);
elseif LeftRightForeHind=='LH'
    if Ssym<0
        symAxisTheta=atan(Ssym);
    else
        symAxisTheta=mod(atan(Ssym),pi)-pi;
    end
    uSlopeLTheta=atan(Sslo);
    Ls2cenTheta= linspace(symAxisTheta,uSlopeLTheta,nSection);
else
    if Ssym<0
        symAxisTheta=atan(Ssym);
    else
        symAxisTheta=mod(atan(Ssym),pi)-pi;
    end
    uSlopeLTheta=atan(Sslo)-pi;
    Ls2cenTheta= linspace(symAxisTheta,uSlopeLTheta,nSection);
end

Ls2cen=-tan(Ls2cenTheta);

end