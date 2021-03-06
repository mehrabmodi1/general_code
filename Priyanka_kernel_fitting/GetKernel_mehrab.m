function [kernelout] = GetKernel_mehrab(starter_kernel, dff_trace, PID_trace)
    
%% start clock
tic

%% 1 : set up the error minimization
Eval_max = 1e+6; 
Iter_max = 1e+6; 
model_fit = @conv_out;  %function conv_out defined below
lb = []; 
ub = [];
options = optimset('MaxFunEvals',Eval_max,'MaxIter',Iter_max);
[kernelout] = lsqcurvefit(model_fit,starter_kernel,PID_trace,dff_trace,lb,ub,options);

%% 2 : define the convolution function (model_fit)
function [dff_trace_predic] = conv_out(starter_kernel,PID_trace)
    %starter_kernel_1 = sgolayfilt(starter_kernel,4,11);
    starter_kernel_1 = starter_kernel;
    x_data = PID_trace(1:end-1, 1);
    %bin_cut(:,1) = PID_trace(end,:);
    z_data = [];
    for i = 1:size(x_data,2)
        M = [];
        M = conv(starter_kernel_1,x_data(:,i));
        M(1,:) = [];
        z_data(1:length(M),i) = M; 
    end
    dff_trace_predic = z_data;
    %dff_trace_predic(dff_trace_predic<0) = 0;
    
    dff_trace_predic = dff_trace_predic(1:size(PID_trace, 1), :);
    
end
dff_trace_predic = [];

%% stop clock
toc

end