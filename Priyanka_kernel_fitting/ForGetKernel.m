function [fitted] = ForGetKernel(temp,mode)
%% 1 : assemble Xdata (odor waveforms), Ydata (FR), Zdata (Predicted FR)

xdata = (temp.odorwfm)';
ydata = (temp.(mode).mean)'; % mode = 'raw' or 'smooth'
sddata = (temp.(mode).sd)';
semdata = (temp.(mode).sem)';
ndata = (temp.(mode).noise)';
blair = temp.baseline;

binmax = max(temp.odorlength);
ydata(binmax+1:end,:) = [];
xdata(binmax+1:end,:) = [];
sddata(binmax+1:end,:) = [];
semdata(binmax+1:end,:) = [];
ndata(binmax+1:end,:) = [];
    
for j = 1:size(ydata,2) 
    % use length of odor waveform to zero off extra bins
    % add last extra bin to keep track of bin counts for each protocol
    bin = temp.odorlength(j,1);
    
    ydata(bin+1:end,j) = 0;
    sddata(bin+1:end,j) = 0;
    semdata(bin+1:end,j) = 0;
    ndata(bin+1:end,j) = 0;    
    
    ydata(binmax+1,j) = bin;
    xdata(binmax+1,j) = bin;
    sddata(binmax+1,j) = bin;
    semdata(binmax+1,j) = bin;
    ndata(binmax+1,j) = bin;    
end

fitted.xdata = xdata;
fitted.ydata = ydata;
fitted.sddata = sddata;
fitted.semdata = semdata;
fitted.ndata = ndata;
fitted.baseline = blair;

end