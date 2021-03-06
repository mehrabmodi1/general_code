function [kernelout] = GetKernel(start_vec,blair,xdata,ydata)
    
%% start clock
tic

%% cautionary step to make sure xdata,ydata and zdata will all be same size
r = max(xdata(end,:)); % largest response length
xdata(r+1:end-1,:) = [];
ydata(r+1:end-1,:) = [];

%% 1 : set up the error minimization
Eval_max = 1e+6; Iter_max = 1e+6; model_fit = @conv_out; lb = []; ub = [];
options = optimset('MaxFunEvals',Eval_max,'MaxIter',Iter_max);
[kernelout] = lsqcurvefit(model_fit,start_vec,xdata,ydata,lb,ub,options);

%% 2 : define the convolution function (model_fit)
function [zdata] = conv_out(start_vec,xdata)
    %start_vec_1 = sgolayfilt(start_vec,4,11);
    start_vec_1 = start_vec;
    x_data = xdata(1:end-1,:);
    bin_cut(:,1) = xdata(end,:);
    z_data = [];
    for i = 1:size(x_data,2)
        M = [];
        M = conv(start_vec_1,x_data(:,i));
        M(1,:) = [];
        M = M + blair;
        z_data(1:length(M),i) = M; % convolved output for single pulse data
    end
    zdata = z_data;
    zdata(zdata<0) = 0;
    for k = 1:size(zdata,2)
        zdata(bin_cut(k,1)+1:end,k) = 0;
    end
    zdata(max(bin_cut)+1:end,:) = [];    
    zdata(end+1,:) = bin_cut;
end
zdata = [];

%% stop clock
toc

end