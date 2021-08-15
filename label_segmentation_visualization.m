function label_matrix_vis=label_segmentation_visualization(label_matrix,color_idx)
        %Create visuallization of segemntation
        label_matrix_visR=zeros(size(label_matrix));
        label_matrix_visG=zeros(size(label_matrix));
        label_matrix_visB=zeros(size(label_matrix));
        
        for idx=1:length(color_idx)
            label_matrix_visR(label_matrix==color_idx(idx,1))=color_idx(idx,2);
            label_matrix_visG(label_matrix==color_idx(idx,1))=color_idx(idx,3);
            label_matrix_visB(label_matrix==color_idx(idx,1))=color_idx(idx,4);     
        end        
        label_matrix_vis=cat(3,label_matrix_visR/255,label_matrix_visG/255,label_matrix_visB/255);
end