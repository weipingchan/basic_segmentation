function [hind_fore_sep,resize_ratio]=get_original_line(mask,json,inline)

    %resize_ratio=json.size(1)/json.resize(1);
    resize_ratio=((size(mask,1)+50)/json.size(1))*(json.size(1)/json.resize(1)); %to overcome the different img size between jpg and matrix; +50 here is for the size of scale bar
    
    l_fh_line=round(inline*resize_ratio);
    label_mask_l=zeros(size(mask,1)-50, size(mask,2)); %Remove the scale bar region
    for pix=1:size(l_fh_line,1)
        label_mask_l(l_fh_line(pix,2),l_fh_line(pix,1))=1;
    end
    label_mask_l_ext=imdilate(label_mask_l,strel('disk',5));
    hind_fore_sep0 = bwskel(logical(label_mask_l_ext),'MinBranchLength',10);
    hind_fore_sep=imdilate(hind_fore_sep0,strel('disk',1));
end