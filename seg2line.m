function infLine=seg2line(LHcropLineSeg,extendLength)
LHcropVector=(LHcropLineSeg(1,:)-LHcropLineSeg(2,:))*extendLength;
infLine=[LHcropLineSeg(1,:)-LHcropVector ; LHcropLineSeg(1,:)+LHcropVector];
end