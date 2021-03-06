%script to make an aligned stack from an analysis results folder.

direc = 'E:\Data\Analysed_data\Manual_ROIs\20190816\60um';
lags = load([direc, '\detailed_xy_lags.mat']);
lags = lags.detailed_lags.frame_lags;
stack = ScanImageTiffReader([direc, '\file_z60_00001.tif']).data();

n_frames = size(stack, 3);
stack_reg = zeros(size(stack, 1), size(stack, 2), size(stack, 3)) + nan;

for frame_n = 1:n_frames
    stack_reg(:, :, frame_n) = translate_stack (squeeze(stack(:, :, frame_n)), [lags(frame_n, 2); lags(frame_n, 1)], nan);
    if rem(frame_n, 50) == 0
        disp(frame_n);
    else
    end
end

PMT_stack = stack_reg(10:(end - 10), 10:(end - 10), 1:2:end);
MPPC_stack = stack_reg(10:(end - 10), 10:(end - 10), 2:2:end);

PMT_im = sqrt(mean(PMT_stack, 3));
MPPC_im = sqrt(mean(MPPC_stack, 3));
%save([direc, '\reg_stack.mat'], 'stack_reg');

save([direc, '\PMT_im.mat'], 'PMT_im');
save([direc, '\PMT_im.mat'], 'MPPC_im');

min_val = min(min(PMT_im));
max_val = max(max(PMT_im));
figure(1)
imagesc(PMT_im, [min_val, (max_val - min_val).*0.8])
colormap('gray')
pbaspect([1 1 1])
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])

min_val = min(min(MPPC_im));
max_val = max(max(MPPC_im));
figure(2)
imagesc(MPPC_im, [min_val, (max_val - min_val).*0.8])
colormap('gray')
pbaspect([1 1 1])
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])

PMT_im = mat2gray(PMT_im);
MPPC_im = mat2gray(MPPC_im);

imwrite(im2uint16(PMT_im), ['E:\Data\Analysed_data\Manual_ROIs\20190816\60um_PMT_aligned_ave.tif']);
imwrite(im2uint16(MPPC_im), ['E:\Data\Analysed_data\Manual_ROIs\20190816\60um_MPPC_aligned_ave.tif']);