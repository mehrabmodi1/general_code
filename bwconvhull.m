function out = bwconvhull(in)
%BWCONVHULL Convex hull of binary image
% bw2 = bwconvhull(bw) computes the convex hull of the input binary image, bw,
% and returns it as another binary image, bw2. bw2 has the same size as bw.
%
% bw must be two-dimensional.

% Steven L. Eddins
% Copyright 2009 The MathWorks, Inc.

if ~islogical(in)
    in = in ~= 0;
end

% Convert input from logical to uint8. This will cause regionprops to treat
% it as a label matrix, and so every foreground pixel will be treated as part of
% the same object.
in = uint8(in);

s = regionprops(in, 'BoundingBox', 'ConvexImage');

% regionprops returns an image only big enough to just cover the bounding box of
% the "object." Compute the row and column indices corresponding to that
% bounding box.
m = s.BoundingBox(4);
n = s.BoundingBox(3);
r1 = s.BoundingBox(2) + 0.5;
c1 = s.BoundingBox(1) + 0.5;
r = (1:m) + r1 - 1;
c = (1:n) + c1 - 1;

out = false(size(in));
out(r,c) = s.ConvexImage;
