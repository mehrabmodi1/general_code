function [stack_al, c_vec, lags] = align_stack(stack_orig, weight_on)
%This function uses matlab's normxcorr2 function to align the frames in 
%stack_orig to the first frame. It also weighs the corrcoef matrix to favour smaller lags when weight_on == 1. 


%test stack
% clear all
% close all
% direc = 'D:\Data\CSHL\test_movt\cell_zoom\cell_zoomed.mat';
% test_stack = load(direc);
% stack_orig = test_stack.ROI_all_pix_mat;


stack_al = zeros(size(stack_orig, 1), size(stack_orig, 2), size(stack_orig, 3)) + nan;
template = squeeze(stack_orig(:, :, 1));
stack_al(:, :, 1) = template;
c_vec = zeros(1, size(stack_orig, 3)) + nan;       %vector of each frame's max corrcoef
c_vec(1) = 1;

%Creating a matrix of size(xcorrcoef matrix) that contains the distance of each
%pixel from the matrix's center. This will be used to weigh the xcorrcoef
%matrix by distance from center (ie lower lags will be favoured).
del = 2.*size(stack_orig, 1) - 1;   %side length of xcorrcoef mat is twice that of stack orig, - 1 for overlapping pix 
dist_mat = zeros(del, del) + nan;
center = [ceil(del./2), ceil(del./2)];
for row_n = 1:del
    for col_n = 1:del
        dist = pdist([center; [row_n, col_n]]);
        dist_mat(row_n, col_n) = dist;
                
    end
    
end
dist_mat(center, center) = 1;
weight_mat = dist_mat./max(max(dist_mat));
weight_mat = abs(weight_mat - 1);

%if weight mat ~= 1, no weighing is done, so weight mat filled with 1's.
if weight_on ~= 1
    weight_mat = zeros(size(weight_mat, 1), size(weight_mat, 2)) + 1;
else
end

%w_powers = [2, 1, 0.5, 0.25];      %range of weight_mat powers tried
w_power = 2.5;                        %this power worked best when assessed manually
curr_weight_mat = weight_mat.^w_power;

%cc_mats_saved = zeros(size(weight_mat, 1), size(weight_mat, 2), 4);

stack_al = zeros(size(stack_orig, 1), size(stack_orig, 2), size(stack_orig, 3)) + nan;
lags = zeros(size(stack_orig, 3), 2) + nan;
for frame_n = 2:size(stack_orig, 3)
    
    
    %for weight_mat_n = 1:4         %used this line of code to try out different weigt_mat powers, 2 worked best
        
        curr_frame = squeeze(stack_orig(:, :, frame_n));

        %using normxcorr2 to calculate offsets
        cc = normxcorr2(template, curr_frame);
        cc = cc.*curr_weight_mat;
        
        %cc_mats_saved(:, :, weight_mat_n) = cc;

        [max_cc, imax] = max(abs(cc(:)));
        c_vec(frame_n) = max_cc;
        
        [ypeak, xpeak] = ind2sub(size(cc),imax(1));
        curr_offset = [ (ypeak-size(template,1)) (xpeak-size(template,2)) ];
        lags(frame_n, :) = curr_offset;
        
        curr_offset(1) = curr_offset(1) .* (- 1);       %fixing sign convention
        

        if curr_offset(1) == 0 && curr_offset(2) == 0
            stack_al(:, :, frame_n) = curr_frame;

        else
            %inserting offset corrected curr frame into stack_al 
            if sign(curr_offset(1)) >= 0 && sign(curr_offset(2)) >= 0
                stack_al( (1):(size(curr_frame, 1) - curr_offset(1) ), (1):(size(curr_frame, 2) - curr_offset(2) ), frame_n)...
                    = curr_frame( (1+ curr_offset(1) ):(size(curr_frame, 1) ), (1+ curr_offset(2) ):(size(curr_frame, 2) ) );    

            elseif sign(curr_offset(1)) < 0 && sign(curr_offset(2)) >= 0
                stack_al(1 + abs(curr_offset(1) ):(size(curr_frame, 1) ), (1):(size(curr_frame, 2) - curr_offset(2) ), frame_n)...
                    = curr_frame( (1):(size(curr_frame, 1) - abs(curr_offset(1)) ), (1 + curr_offset(2) ):(size(curr_frame, 2) ) );

            elseif sign(curr_offset(1)) >= 0 && sign(curr_offset(2)) < 0
                stack_al( (1):(size(curr_frame, 1) - curr_offset(1) ), (1 + abs(curr_offset(2)) ):(size(curr_frame, 2) ), frame_n)...
                    = curr_frame(1 + (curr_offset(1) ):(size(curr_frame, 1) ), (1):(size(curr_frame, 2) - abs(curr_offset(2) ) ) );

            elseif sign(curr_offset(1)) < 0 && sign(curr_offset(2)) < 0
                stack_al(1 + abs(curr_offset(1) ):(size(curr_frame, 1) ), (1 + abs(curr_offset(2) ) ):(size(curr_frame, 2) ), frame_n)...
                    = curr_frame( (1):(size(curr_frame, 1) - abs(curr_offset(1) ) ), (1):(size(curr_frame, 2) - abs(curr_offset(2) ) ) );


            else
            end

            
            
        end
%         figure(1)
%         subplot(2, 2, 1)
%         imagesc(template)
% 
%         subplot(2, 2, 2)
%         imagesc(curr_frame)
%         title(['weight mat power ' int2str(w_powers(weight_mat_n)) ])
%         
%         subplot(2, 2, 3)
%         imagesc(stack_al(:, :, frame_n))
%         title(int2str(curr_offset))
%         
%         subplot(2, 2, 4)
%         imagesc(cc)
%         drawnow
% 
%         del = input('press enter');
%         stack_al(:, :, frame_n) = nan;
end
%        figure(2)
%         subplot(2, 2, 1)
%         imagesc(cc_mats_saved(:, :, 1))
%         title(int2str(curr_offsets_saved(1, :)))
%         
%         subplot(2, 2, 2)
%         imagesc(cc_mats_saved(:, :, 2))
%         title(int2str(curr_offsets_saved(2, :)))
% 
%         subplot(2, 2, 3)
%         imagesc(cc_mats_saved(:, :, 3))
%         title(int2str(curr_offsets_saved(3, :)))
% 
%         subplot(2, 2, 4)
%         imagesc(cc_mats_saved(:, :, 4))
%         title(int2str(curr_offsets_saved(4, :)))
%     
%         
%     
end
