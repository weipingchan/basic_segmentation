function [NewImg,NewOrigion]=rmBlankRegion(inImg,bufferWidth)
    bwInImg=imbinarize(grayImg(inImg));
    BbwIn=regionprops(uint8(bwInImg),'BoundingBox'); %The center of the bounding box
    boundingBoxbwIn=round(BbwIn.BoundingBox);
    %bufferWidth=100;
    if boundingBoxbwIn(2)-bufferWidth<1 upbond=1; else upbond=boundingBoxbwIn(2)-bufferWidth;, end;
    if boundingBoxbwIn(2)+boundingBoxbwIn(4)+bufferWidth>size(inImg,1) lowbond=size(inImg,1);, else lowbond=boundingBoxbwIn(2)+boundingBoxbwIn(4)+bufferWidth;, end;
    if boundingBoxbwIn(1)-bufferWidth<1 leftbond=1;, else leftbond=boundingBoxbwIn(1)-bufferWidth;, end;
    if boundingBoxbwIn(1)+boundingBoxbwIn(3)+bufferWidth>size(inImg,2) rightbond=size(inImg,2);, else rightbond=boundingBoxbwIn(1)+boundingBoxbwIn(3)+bufferWidth;, end;
    NewImg = inImg(upbond:lowbond,leftbond:rightbond,:);
    NewOrigion=[leftbond,upbond];
end