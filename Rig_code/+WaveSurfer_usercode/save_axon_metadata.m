function metadata = save_axon_metadata()
%This function uses MulticlampTelegraph to get amp electrode parameters
%from the amplifier and adds them to the structure that will be saved by
%WaveSurfer in it's data file.
%Mehrab Modi 20180207

%getting electrode id
a = ws.dabs.axon.MulticlampTelegraph('getAllElectrodeIDs');

%getting electrode metadata - assumes there's only one electrode
inf = ws.dabs.axon.MulticlampTelegraph('getElectrodeState', a);

keyboard