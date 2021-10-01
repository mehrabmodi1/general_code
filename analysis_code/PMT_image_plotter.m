clear all
close all

%tif_path = 'C:\Data\Data\Analysed_data\Suite2P_results\20210615\fly1_c739_opGC6f\1\file_00001.tif';                 %KCs

tif_path = 'C:\Data\Data\Analysed_data\Manual_ROI_results\20200805\fly2_MBONG2Ap1_simp_pairing_Berry\odor_trs_00001.tif';      %MBON G2A'1
%tif_path = 'C:\Data\Data\Analysed_data\Manual_ROI_results\20200805\fly3_MBONG2Ap1_simp_pairing_Berry\odor_trs_00002.tif';      %MBON G2A'1

script_name = mfilename;

if isempty(strcmp(tif_path, 'c739')) == 0
    pix_per_mic = 0.22;     %pixels per micron, KC datasets have zoom = 9
else
    pix_per_mic = 0.18;     %pixels per micron, MBON datasets have zoom = 11
end


stim_frs = [101, 151];
n_cells = 2;
frame_time = 0.091;

stack_obj = ScanImageTiffReader(tif_path);
[frame_time, zoom, n_chans, PMT_offsets, lpower] = SI_tif_info(stack_obj);
stack = stack_obj.data();
stack = permute(stack,[2 1 3]);
stack = double(stack);
if n_chans == 2
    stack = stack(:, :, 1:2:end);
else
end

PMT_stack = align_image_rows(stack, []);        %accounting for pixel row offsets due to imperfect turn-around time compenstaion

%manually measured background values (PMT offset)
bk_val_PMT = 110;
%subtracting background
PMT_stack = PMT_stack - bk_val_PMT;

%mean frame
ave_im_PMT = mean(PMT_stack, 3);

%baseline frames
base_fr_PMT = mean(PMT_stack(:, :, 1:(stim_frs(1) - 10)), 3, 'omitnan');
base_fr_PMT_bar = add_scale_bar(base_fr_PMT, 0.18, 10, 1);
base_fr_PMT(base_fr_PMT < 70) = 0;                              %thresholding away bacground to hide it in dF/F image

figure(1)
PMT_im = real(base_fr_PMT_bar.^0.5);
imagesc(PMT_im, [0, max(max(PMT_im)).*0.8])
colormap('gray')
pbaspect([1 1 1])
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])


%resp period frames
resp_fr_PMT = mean(PMT_stack(:, :, stim_frs(1):stim_frs(2)), 3, 'omitnan');

%dF/F frames
dFF_fr_PMT = (resp_fr_PMT - base_fr_PMT)./base_fr_PMT;
dFF_fr_PMT(isinf(dFF_fr_PMT)) = 0;

figure(4)
imagesc(dFF_fr_PMT, [0, 2.5])
colormap('gray')
pbaspect([1 1 1])
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])
%colorbar
pbaspect([1 1 1])

[diff_fr, sig_resp_frame, sig_resp_frame_n, baseline_fr, stack] = highlight_resp_pix(3, PMT_stack, stim_frs, 1, 0.099, pix_per_mic, 0);

baseline_stack = repmat(baseline_fr, 1, 1, size(stack, 3));


%extracting response traces for manually chosen cell ROIs
ROI_mat = [];
trace_mat = zeros(size(stack, 3), n_cells) + nan;
for cell_n = 1:n_cells
    figure(3)
    BW = roipoly;
    ROI_mat = cat(3, ROI_mat, BW);
    curr_pix = find(BW == 1);
    
    for frame_n = 1:size(stack, 3)
        curr_frame = squeeze(stack(:, :, frame_n));
        trace_mat(frame_n, cell_n) = mean(mean(curr_frame(curr_pix), 'omitnan'), 'omitnan');
    end
end

%computing dF/F values for traces
F0 = mean(trace_mat(10:(stim_frs(1) - 10), :), 'omitnan');
F0 = repmat(F0, size(trace_mat, 1), 1);
trace_mat = (trace_mat - F0)./F0;

%plotting extracted traces
cascade_plot(5, trace_mat, [0.6, 0.6, 0.6], 1, 0, 4, 4, 5, 0, 1);
set_xlabels_time(5, frame_time, 2);
add_stim_bar(5, stim_frs, [0.65, 0.65, 0.65] )

