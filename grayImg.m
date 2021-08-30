function gray=grayImg(inimg)
%Convert RGB images into gray scale; the gray image stays the same
    [~, ~, chab]=size(inimg);
    if chab>1
        gray=rgb2gray(inimg);
    else
        gray=inimg;
    end
end