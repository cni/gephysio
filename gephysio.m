function gephysio(phys_dir, out_dir, scan_TR, nvolumes, plot_flag, plot_start, plot_window, ppg_dt, resp_dt, slice_timing)
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
%   ppg_dt              PPG data sampling window (ms)
%   resp_dt             respiratory data sampling window (ms)
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
    ppg_dt      = str2double(ppg_dt);
    resp_dt     = str2double(resp_dt);
end
if ~exist('ppg_dt', 'var') || isempty(ppg_dt), ppg_dt = 10; end
if ~exist('resp_dt', 'var') || isempty(resp_dt), resp_dt = 40; end
 
param.PPG.dt  = ppg_dt;        % PPG data samples at 10ms
param.RESP.dt = resp_dt;       % Respiratory data samples  at 40 ms
dirlist = dir(fullfile(phys_dir, 'PPGData*'));
param.PPG.wave.fn  = fullfile(dirlist.folder, dirlist.name);
dirlist = dir(fullfile(phys_dir, 'PPGTrig*'));
param.PPG.trig.fn  = fullfile(dirlist.folder, dirlist.name);
dirlist = dir(fullfile(phys_dir, 'RESPData*'));
param.RESP.wave.fn = fullfile(dirlist.folder, dirlist.name);
dirlist = dir(fullfile(phys_dir, 'RESPTrig*'));
param.RESP.trig.fn = fullfile(dirlist.folder, dirlist.name);

param.TR             = scan_TR;         % milliseconds
param.nvols          = nvolumes;
param.scan_duration  = param.TR * param.nvols;   
if exist('slice_timing', 'var') && ~isempty(slice_timing)
    param.sliceTiming    = slice_timing;
    [~,param.sliceOrder] = sort(param.sliceTiming);
end

% Read the physio data and find peaks

physioType = {'PPG', 'RESP'};

% First read in the waveform data to determine the total length of the recording,
for p = 1:numel(physioType)
    [param.(physioType{p}).wave.data_raw, param.(physioType{p}).wave.data_sync, ...
        param.(physioType{p}).wave.t_raw, param.(physioType{p}).wave.t_sync, param.scanStart] = ...
        physioRead(param.(physioType{p}).wave.fn, param.(physioType{p}).dt, param.scan_duration, 'wave');    
    % scale the sync waveform to [0, 1]
    param.(physioType{p}).wave.data_sync = double(param.(physioType{p}).wave.data_sync - min(param.(physioType{p}).wave.data_sync)) / ...
        (max(param.(physioType{p}).wave.data_sync) - min(param.(physioType{p}).wave.data_sync));

end
% then trim the trigger data accordingly
for p = 1:numel(physioType)
    [param.(physioType{p}).trig.data_raw, param.(physioType{p}).trig.data_sync,~,~,~] = ...
        physioRead(param.(physioType{p}).trig.fn, param.(physioType{p}).dt, param.scan_duration, 'trig', param.scanStart);
end

% Filter signal and find peaks in PPG and RESP data using Dora's function
for p = 1:numel(physioType)
    [param.(physioType{p}).trig.data_filt, param.(physioType{p}).wave.data_filt] = ...
        physioPeaks(param.(physioType{p}).wave.data_sync, param.(physioType{p}).dt);
end
% set up a flag array for the the detected peaks
for p = 1:numel(physioType)
    param.(physioType{p}).trig.filter_flag = zeros(size(param.(physioType{p}).wave.t_sync));
    [~,loc] = ismember(round(param.(physioType{p}).trig.data_filt), param.(physioType{p}).wave.t_sync);
    param.(physioType{p}).trig.filter_flag(loc) = 1;
end
    
%% Plot a segment of the data
if exist('plot_flag', 'var') && plot_flag

if ~exist('plot_window', 'var') || isempty(plot_window), plot_window = 10000; end   % Duration of physio data to plot for inspection (10 seconds)
if ~exist('plot_start', 'var') || isempty(plot_start),   plot_start  = 0;  end
if plot_start<param.RESP.dt, plot_start=param.RESP.dt; end
if plot_start>param.scan_duration-plot_window, plot_start=param.scan_duration-plot_window; end

for p = 1:numel(physioType)
    h = figure(p*100);
    set(h, 'Visible', 'off');
    subplot(3,1,1);
    plot_start = 0; plot_end = plot_start + plot_window;
    r = (plot_start < param.(physioType{p}).wave.t_sync) & (param.(physioType{p}).wave.t_sync < plot_end);
    plot(param.(physioType{p}).wave.t_sync(r), param.(physioType{p}).wave.data_sync(r),'b'); hold on;
    plot(param.(physioType{p}).wave.t_sync(r), param.(physioType{p}).wave.data_filt(r),'r--', 'LineWidth',1); 
    
    trig_in_window = param.(physioType{p}).trig.data_filt(...
        (plot_start < param.(physioType{p}).trig.data_filt) & ...
        (param.(physioType{p}).trig.data_filt < plot_end));
    plot(trig_in_window, ones(size(trig_in_window))*max(param.(physioType{p}).wave.data_sync),'ro');
    title(sprintf('Beginning of scan: %d - %d s', floor(plot_start/1000), floor(plot_end/1000)));
    xlabel('time (ms)'); ylabel(sprintf('%s signal', physioType{p})); xlim([plot_start, plot_end]);
    legend({'recorded signal','filtered signal','trigger'}, 'Location','south','NumColumns',3,'FontSize',10);legend('boxoff'); 
    
    subplot(3,1,2);
    plot_start = param.scan_duration / 2; plot_end = plot_start + plot_window;
    r = (plot_start < param.(physioType{p}).wave.t_sync) & (param.(physioType{p}).wave.t_sync < plot_end);
    plot(param.(physioType{p}).wave.t_sync(r), param.(physioType{p}).wave.data_sync(r),'b'); hold on;
    plot(param.(physioType{p}).wave.t_sync(r), param.(physioType{p}).wave.data_filt(r),'r--', 'LineWidth',1); 
    
    trig_in_window = param.(physioType{p}).trig.data_filt(...
        (plot_start < param.(physioType{p}).trig.data_filt) & ...
        (param.(physioType{p}).trig.data_filt < plot_end));
    plot(trig_in_window, ones(size(trig_in_window))*max(param.(physioType{p}).wave.data_sync),'ro');
    title(sprintf('Middle of scan: %d - %d s', floor(plot_start/1000), floor(plot_end/1000)));
    xlabel('time (ms)'); ylabel(sprintf('%s signal', physioType{p})); xlim([plot_start, plot_end]);
    legend({'recorded signal','filtered signal','trigger'}, 'Location','south','NumColumns',3,'FontSize',10);legend('boxoff'); 

    subplot(3,1,3);
    plot_start = param.scan_duration - plot_window; plot_end = param.scan_duration;
    r = (plot_start < param.(physioType{p}).wave.t_sync) & (param.(physioType{p}).wave.t_sync < plot_end);
    plot(param.(physioType{p}).wave.t_sync(r), param.(physioType{p}).wave.data_sync(r),'b'); hold on;
    plot(param.(physioType{p}).wave.t_sync(r), param.(physioType{p}).wave.data_filt(r),'r--', 'LineWidth',1); 
    
    trig_in_window = param.(physioType{p}).trig.data_filt(...
        (plot_start < param.(physioType{p}).trig.data_filt) & ...
        (param.(physioType{p}).trig.data_filt < plot_end));
    plot(trig_in_window, ones(size(trig_in_window))*max(param.(physioType{p}).wave.data_sync),'ro');
    title(sprintf('End of scan: %d - %d s', floor(plot_start/1000), floor(plot_end/1000)));
    xlabel('time (ms)'); ylabel(sprintf('%s signal', physioType{p})); xlim([plot_start, plot_end]);
    legend({'recorded signal','filtered signal','trigger'}, 'Location','south','NumColumns',3,'FontSize',10);legend('boxoff'); 

    % Programatically move the Legend
%     hL = legend({'recorded signal','filtered signal','trigger'},'NumColumns',2);
%     newPosition = [0.7 0.9 0.2 0.2];
%     newUnits = 'normalized';
%     set(hL,'Position', newPosition,'Units', newUnits, 'Box', 'off');
    
    fname = fullfile(out_dir, sprintf('%s_SampSig', physioType{p}));
    print(fname, '-dpdf', '-fillpage');

end

close all;  % close all figures
end

%% Write out simple summary of the data
% average heart beat / respiration and standard deviation in millisecond
for p = 1:numel(physioType)
    param.(physioType{p}).trig.dtrig = diff(param.(physioType{p}).trig.data_filt);
    param.(physioType{p}).trig.std = std(param.(physioType{p}).trig.dtrig); 
    param.(physioType{p}).trig.mean = mean(param.(physioType{p}).trig.dtrig(abs(param.(physioType{p}).trig.dtrig - ...
        mean(param.(physioType{p}).trig.dtrig)) < 3 * std(param.(physioType{p}).trig.dtrig))); 
end

% Write filted signal and detected peaks to text files. Save the filtered
% signal as integers and peaks are rounded to millisecond
for p = 1:numel(physioType)
    dlmwrite(fullfile(out_dir, sprintf('%s_FltData.csv',physioType{p})), [param.(physioType{p}).wave.t_sync, param.(physioType{p}).wave.data_filt, param.(physioType{p}).trig.filter_flag], 'precision', '%d');
    dlmwrite(fullfile(out_dir, sprintf('%s_FltTrig.csv',physioType{p})), int64(param.(physioType{p}).trig.data_filt), 'precision', '%d');
    fid = fopen(fullfile(out_dir, sprintf('%s_Stat.csv',physioType{p})), 'w');
    fprintf(fid, 'Scan duration (s): %.1f\n', param.scan_duration/1000);
    fprintf(fid, 'Number of peaks detected: %d\n', length(param.(physioType{p}).trig.data_filt));
    fprintf(fid, 'Peak interval standard deviation (s): %.3f, mean (s): %.3f\n', param.(physioType{p}).trig.std/1000, param.(physioType{p}).trig.mean/1000);
    fprintf(fid, 'Peak intervals (s):\n');
    for i = 1:length(param.(physioType{p}).trig.dtrig)
        fprintf(fid, '%.3f\n', param.(physioType{p}).trig.dtrig(i)/1000);
    end
    fclose(fid);
end


% if verLessThan('matlab', '9.10')
%     str = jsonencode(param);
% else
%     str = jsonencode(param, 'PrettyPrint', true);
% end
% fid = fopen(fullfile(out_dir, 'physio.json'), 'w');
% fprintf(fid, str);
% fclose(fid);

end
