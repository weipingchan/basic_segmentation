function infLine=seg2line(LHcropLineSeg,extendLength)
%Convert a segment into a line vector
LHcropVector=(LHcropLineSeg(1,:)-LHcropLineSeg(2,:))*extendLength;
infLine=[LHcropLineSeg(1,:)-LHcropVector ; LHcropLineSeg(1,:)+LHcropVector];
end