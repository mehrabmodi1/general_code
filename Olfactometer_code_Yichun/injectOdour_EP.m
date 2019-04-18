function injectOdour_EP(vial)
if nargin<2, duration=[]; end

AIR_VIAL=13; %This is fixed unless we re-build the olfactometer

if vial==AIR_VIAL
    return
end

 FlipValve_EP({sprintf('Vial%d',vial),'NO'},0)
