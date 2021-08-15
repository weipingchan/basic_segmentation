function beltWpar=findbeltWpar(mask,forehindCornerL,forehindCornerR,realCen,boundingBox)
try
    areaSum=sum(mask,1);
    axis0=[1:size(mask,2)]-round(realCen(1));
    ferqSeries0=[axis0;smoothdata(areaSum,'sgolay')].';
    areaBoundaryLengthR=findbeltWpar0(ferqSeries0,'R',forehindCornerR,realCen);
    areaBoundaryLengthL=findbeltWpar0(ferqSeries0,'L',forehindCornerL,realCen);
    beltWpar0=ceil((areaBoundaryLengthR+areaBoundaryLengthL+20)*2/boundingBox(3)*100)/100;
catch
    beltWpar0=0.2;
end
%make sure the minimal beltWpar greater than 0.2
if beltWpar0<0.2 || beltWpar0>0.4
    beltWpar =0.2;
else
    beltWpar =beltWpar0;
end
end