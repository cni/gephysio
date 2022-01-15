function gephysio(phys_dir, out_dir, scan_TR, nvolumes, plot_flag, plot_start, plot_window, slice_timing)
%
% Parse GE physio recordings, smooth the pulse oximetry and respiratory
% waveform, and find peaks in the waveform.
%
% Input
%   phys_dir            directory containing the physio files
%   out_dir             output directory
%   scan_TR             scan TR, in millisecond
%   nvolumes            number of volumes, integer
%   plot_flag           flag for plotting a segment of the physio data
%   plot_start
%   plot_window
%   slice_timing        (not used)

if ~exist('out_dir','var') || isempty(out_dir)
    out_dir = phys_dir;
end
% Convert from string to number if running as compiled Matlab binary
if isdeployed 
    scan_TR     = str2double(scan_TR);
    nvolumes    = str2double(nvolumes);
    plot_flag   = str2double(plot_flag);
    plot_start  = str2double(plot_start);
    plot_window = str2double(plot_window);
end

param.ppg.dt  = 10;       % PPG data samples at 10ms
param.resp.dt = 40;       % Respiratory data samples  at 40 ms
dirlist = dir(fullfile(phys_dir, 'PPGData*'));
param.ppg.wave.fn  = fullfile(dirlist.folder, dirlist.name);
dirlist = dir(fullfile(phys_dir, 'PPGTrig*'));
param.ppg.trig.fn  = fullfile(dirlist.folder, dirlist.name);
dirlist = dir(fullfile(phys_dir, 'RESPData*'));
param.resp.wave.fn = fullfile(dirlist.folder, dirlist.name);
dirlist = dir(fullfile(phys_dir, 'RESPTrig*'));
param.resp.trig.fn = fullfile(dirlist.folder, dirlist.name);

param.TR             = scan_TR;         % milliseconds
param.nvols          = nvolumes;
param.scan_duration  = param.TR * param.nvols;   
if exist('slice_timing', 'var') && ~isempty(slice_timing)
    param.sliceTiming    = slice_timing;
    [~,param.sliceOrder] = sort(param.sliceTiming);
end

% Read the physio data and find peaks

physioType = {'ppg', 'resp'};

% First read in the waveform data to determine the total length of the recording,
% then trim the trigger data accordingly
for p = 1:numel(physioType)
    [param.(physioType{p}).wave.data_raw, param.(physioType{p}).wave.data_sync, ...
        param.(physioType{p}).wave.t_raw, param.(physioType{p}).wave.t_sync, param.scanStart] = ...
        physioRead(param.(physioType{p}).wave.fn, param.(physioType{p}).dt, param.scan_duration, 'wave');
end

for p = 1:numel(physioType)
    [param.(physioType{p}).trig.data_raw, param.(physioType{p}).trig.data_sync,~,~,~] = ...
        physioRead(param.(physioType{p}).trig.fn, param.(physioType{p}).dt, param.scan_duration, 'trig', param.scanStart);
end

% Filter signal and find peaks in ppg and resp data using Dora's function
for p = 1:numel(physioType)
    [param.(physioType{p}).trig.data_filt, param.(physioType{p}).wave.data_filt] = ...
        physioPeaks(param.(physioType{p}).wave.data_sync, param.(physioType{p}).dt);
end


% Plot a segment of the data
if exist('plot_flag', 'var') && plot_flag
h = figure(101);
set(h, 'Visible', 'off');

if ~exist('plot_window', 'var') || isempty(plot_window), plot_window = 10000; end   % Duration of physio data to plot for inspection (10 seconds)
if ~exist('plot_start', 'var') || isempty(plot_start),   plot_start  = 0;  end
if plot_start<param.resp.dt, plot_start=param.resp.dt; end
if plot_start>param.scan_duration-plot_window, plot_start=param.scan_duration-plot_window; end
plot_end    = plot_start + plot_window;

for p = 1:numel(physioType)
    subplot(2,1,p);
    r = (plot_start < param.(physioType{p}).wave.t_sync) & (param.(physioType{p}).wave.t_sync < plot_end);
    plot(param.(physioType{p}).wave.t_sync(r), param.(physioType{p}).wave.data_sync(r),'b'); hold on;
    plot(param.(physioType{p}).wave.t_sync(r), param.(physioType{p}).wave.data_filt(r),'r:', 'LineWidth',1); 
    
    getrig_in_window = param.(physioType{p}).trig.data_sync(...
        (plot_start < param.(physioType{p}).trig.data_sync) & ...
        (param.(physioType{p}).trig.data_sync < plot_end));
    trig_in_window = param.(physioType{p}).trig.data_filt(...
        (plot_start < param.(physioType{p}).trig.data_filt) & ...
        (param.(physioType{p}).trig.data_filt < plot_end));
    plot(getrig_in_window, ones(size(getrig_in_window))*max(param.(physioType{p}).wave.data_sync),'bx');
    plot(trig_in_window, ones(size(trig_in_window))*max(param.(physioType{p}).wave.data_sync),'ro');
    xlabel('time (ms)'); ylabel(sprintf('%s recorded data', physioType{p})); xlim([plot_start, plot_end]);
    legend({'recorded waveform','smoothed waveform','recorded trigger','smoothed trigger'}, 'Location','south','NumColumns',2);legend('boxoff'); 
end

pdfname = fullfile(out_dir, ['physio_' num2str(plot_start) '-' num2str(plot_end) 'ms.pdf']);
saveas(h, pdfname);
end

% Write out simple summary of the data
% average heart beat / respiration rate in Hz
for p = 1:numel(physioType)
    dtrig = diff(param.(physioType{p}).trig.data_filt);
    param.(physioType{p}).trig.average_rate = 1000 / mean(dtrig(abs(dtrig - mean(dtrig)) < 3 * std(dtrig))); 
end

if verLessThan('matlab', '9.10')
    str = jsonencode(param);
else
    str = jsonencode(param, 'PrettyPrint', true);
end
fid = fopen(fullfile(out_dir, 'physio.json'), 'w');
fprintf(fid, str);
fclose(fid);

end
