function [fwhm x1 x2] = lin_interp_fwhm(y_vals, pki, frac)

%syntax [fwhm x1 x2] = lin_interp_fwhm(yvals xvals pki) This function
%calculates the full width at half max curve defined by yvals and xvals -
%where the pk of the curve is the highest yval (but specified by pki). The
%outputs x1 and x2 are values in the space defined by xvals such that they
%best represent where half y would be half the value at pki, on either side of it as defined by
%linear interpolation. fwhm is the full-width at half maximum - the
%difference between x1 and x2.
% frac is the fraction of the maximum value used to calculate the width of
% the peak - default is 0.5


%test sample - comment out first line and un-comment the next three to test function
% y_vals = [0, 1, 3, 6.5, 7, 3.6, 5.5, 1];
% x_vals = 1:1:length(y_vals);
%plot(y_vals)

if nargin == 1
    [pk pki] = max(y_vals);
else
end

left_yvals = y_vals(1:pki);
right_yvals = y_vals(pki:length(y_vals));

if nargin == 2
    halfmax = y_vals(pki)./2;
elseif nargin == 3
    halfmax = y_vals(pki).*frac;
else
end

%finding points nearest to pki in y_vals where sign changes on subtracting halfmax
diff_vec_l = left_yvals - halfmax;  %left half first (y_vals(i) where i < pki)
signs = sign(diff_vec_l);
signs = diff(signs);
changei = find(signs ~= 0);
changei = max(changei);
lefti = [changei, changei + 1];

diff_vec_r = right_yvals - halfmax; %right half
signs = sign(diff_vec_r);
signs = diff(signs);
changei = find(signs ~= 0);
changei = min(changei);
righti = [changei, changei + 1] + pki - 1;

if isempty(righti) == 1 | isempty(lefti) == 1
    fwhm = nan;
    x1 = nan;
    x2 = nan;
    disp('peak not at pki')
else

    clear signs
    clear changei
    clear diff_vec_r
    clear diff_vec_l
    clear left_yvals
    clear right_yvals


    %fitting lines to the two pairs of points and finding x_val corress to
    %y_val = halfmax on these two lines

    %left of pk
    m = (y_vals(lefti(1, 2)) - y_vals(lefti(1, 1)) )./(lefti(1, 2) - lefti(1, 1));
    c = y_vals(lefti(1, 1)) - lefti(1, 1).*m;

    x1 = (halfmax - c)./m;      %obtained from y = mx + c


    %right of pk
    m = (y_vals(righti(1, 2)) - y_vals(righti(1, 1)) )./(righti(1, 2) - righti(1, 1));
    c = y_vals(righti(1, 1)) - righti(1, 1).*m;

    x2 = (halfmax - c)./m;


    %fwhm calculated as x2 - x1
    fwhm = x2 - x1;
end