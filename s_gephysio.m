%% s_gephysio
%
% Author: Hua Wu
%
% Illustrates how to download data from the iris project, part of the
% LEANEW1 group 
% 
% Then we read the data.
% Then we process the data using the DH methods.
% Then we save the file.
%
% We will convert this into a Gear at some point.
%

% General information so we can find the physio data
group         = 'leanew1';
project_label = 'iris';
subj_label    = 'IRIS021';
session_label = '22776';
acq_label     = 'resting_2';

%% Here are the scitran commands for finding and downloading the data

% see bbResponse2physio in BrainBeats

% Let's try with scitran for a while.
% st = scitran('cni');

% Or use Flywheel MATLAB SDK. fw = flywheel.Flywheel(apikey)
% addpath('~/Documents/MATLAB/'); fw = setup_flywheel;
if exist('fw', 'var') && isa(fw, 'flywheel.Flywheel')
    userInfo = fw.getCurrentUser;
    if ~isempty(userInfo)
        fprintf('Continue with current Flywheel user: %s \n', userInfo.id);
    else
        warning('Empty Flywheel user info. \n');
    end
else
    fprintf('Connect to Flywheel... \n');
    apikey = input('Please enter apikey: ', 's');
    fw = flywheel.Flywheel(apikey);
    userInfo = fw.getCurrentUser;
    fprintf('Logged into Flywheel as %s\n', userInfo.id);
end

datadir = fullfile(bbPath,'local');

%% Download files from Flywheel using scitran

% Using lookup to find data. Examples of lookup methods that include the subject
lookup = fullfile(group,project_label,subj_label,session_label,acq_label);
thisAcq = fw.lookup(lookup);

% See the file names
for f=1:numel(thisAcq.files) 
    % fprintf([thisAcq.files{f}.name '\n']); 
    if contains(thisAcq.files{f}.name, 'gephysio')
        gephysioZip = thisAcq.files{f};
    elseif contains(thisAcq.files{f}.name, 'dicom')
        dcmFile = thisAcq.files{f};
    elseif contains(thisAcq.files{f}.name, 'nii')
        niFile = thisAcq.files{f};
    end
end
% Download and unzip the physio files
gephysioFile = fullfile(datadir, gephysioZip.name);
gephysioZip.download(gephysioFile);
physiofiles = unzip(gephysioFile, datadir);
niFile.download(fullfile(datadir, niFile.name));
rsfmri = niftiread(fullfile(datadir, niFile.name));

% set physio data filenames and sampling rates
param.ppg.dt  = 10;       % PPG data samples at 10ms
param.resp.dt = 40;       % Respiratory data samples  at 40 ms
param.ppg.wave.fn  = physiofiles{contains(physiofiles,'PPGData')};
param.ppg.trig.fn  = physiofiles{contains(physiofiles,'PPGTrig')};
param.resp.wave.fn = physiofiles{contains(physiofiles,'RESPData')};
param.resp.trig.fn = physiofiles{contains(physiofiles,'RESPTrig')};
% Get the acquisition parameters from the dicom, nifti file metadata (info)
param.TR             = dcmFile.info.RepetitionTime;
param.nvols          = dcmFile.info.NumberOfTemporalPositions;
param.scan_duration  = param.TR * param.nvols;   % milliseconds
param.sliceTiming    = niFile.info.SliceTiming;
[~,param.sliceOrder] = sort(param.sliceTiming);

%% Read the physio data and find peaks

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

%% Plot a segment of the data

plot_window = 10000;    % Duration of physio data to plot for inspection (10 seconds)
plot_start  = 4000; 
if plot_start<param.resp.dt, plot_start=param.resp.dt; end
if plot_start>param.scan_duration-plot_window, plot_start=param.scan_duration-plot_window; end
plot_end    = plot_start + plot_window;

figure; 
for p = 1:numel(physioType)
    subplot(2,1,p);
    r = (plot_start < param.(physioType{p}).wave.t_sync) & (param.(physioType{p}).wave.t_sync < plot_end);
    plot(param.(physioType{p}).wave.t_sync(r), param.(physioType{p}).wave.data_sync(r),'b'); hold on;
    plot(param.(physioType{p}).wave.t_sync(r), param.(physioType{p}).wave.data_filt(r),'r:', 'LineWidth',1); 
    
    trig_in_window = param.(physioType{p}).trig.data_filt(...
        (plot_start < param.(physioType{p}).trig.data_filt) & ...
        (param.(physioType{p}).trig.data_filt < plot_end));
    plot(trig_in_window, ones(size(trig_in_window))*max(param.(physioType{p}).wave.data_sync),'ro');
    getrig_in_window = param.(physioType{p}).trig.data_sync(...
        (plot_start < param.(physioType{p}).trig.data_sync) & ...
        (param.(physioType{p}).trig.data_sync < plot_end));
    plot(getrig_in_window, ones(size(getrig_in_window))*max(param.(physioType{p}).wave.data_sync),'bx');
    xlabel('Time (ms)'); ylabel(sprintf('%s recorded data', physioType{p})); xlim([plot_start, plot_end]);
end

%% find ppg trigger for each fmri volume

param.rsfmri.t = (0:param.TR:(param.scan_duration - param.TR));
hbeats = param.ppg.trig.data_filt;
for sl = 1:size(rsfmri, 3)
    for ii = 1:numel(param.rsfmri.t)
        tt = param.rsfmri.t(ii) + param.sliceTiming(sl);    % correct for sliceTiming for each slice
        hbeats_before = (hbeats < tt);
        this_beat = hbeats(find(hbeats_before, 1, 'last'));
        if isempty(this_beat)
            param.rsfmri.beat_tdiff(sl,ii) = NaN;
        else
            param.rsfmri.beat_tdiff(sl,ii) = tt - this_beat;
        end
    end
end

% regress the fmri time series 
dummy_cycles = 3;
slices = 1:size(rsfmri,3);
tpoints = dummy_cycles : size(rsfmri,4);
signal_raw_v = single(reshape(rsfmri(:,:,slices,tpoints),[],numel(tpoints)));
signal_reg_v = zeros(size(signal_raw_v));
for voxel = 1:size(signal_raw_v,1)
    signal = signal_raw_v(voxel,:);
    coeff = polyfit(param.rsfmri.t(tpoints), signal, 2);
    signal_fit = polyval(coeff, param.rsfmri.t(tpoints));
    signal_reg_v(voxel,:) = signal - signal_fit;
end

rsfmri_reg = reshape(signal_reg_v, size(rsfmri(:,:,slices,tpoints)));

sl = 20;
x = 60;
y = 39;
td = param.rsfmri.beat_tdiff(sl,tpoints);
tbin_size = 100; % 100ms 
tbin_shift = 10;
nsteps = tbin_size/tbin_shift - 1;
nbins = ceil(max(td)/tbin_size);
for jj = 1:(tbin_size/tbin_shift - 1)
tbins = (jj-1)*tbin_shift:tbin_size:max(td);
for ii = 1:numel(tbins)
    idx = (tbins(ii) < td) & (td < tbins(ii)+tbin_size);
    signal_mean(jj+(ii-1)*nbins) = mean(rsfmri_reg(x,y,sl,idx));
end
end

figure; plot(param.rsfmri.beat_tdiff(sl,tpoints),squeeze(rsfmri_reg(x,y,sl,:)),'.'); hold on;
plot(1:tbin_shift:numel(signal_mean)*tbin_shift,signal_mean,'ro')


%% find fmri slice closest to each ppg trigger

param.rsfmri.t = (0:param.TR:(param.scan_duration - param.TR));
hbeats = param.ppg.trig.data_filt;
for sl = 1:size(rsfmri, 3)
    for ii = 1:numel(param.rsfmri.t)
        tt = param.rsfmri.t(ii) + param.sliceTiming(sl);    % correct for sliceTiming for each slice
        hbeats_before = (hbeats < tt);
        this_beat = hbeats(find(hbeats_before, 1, 'last'));
        if isempty(this_beat)
            param.rsfmri.beat_tdiff(sl,ii) = NaN;
        else
            param.rsfmri.beat_tdiff(sl,ii) = tt - this_beat;
        end
    end
end



%% Write out simple summary of the data
% average heart beat / respiration rate in Hz
for p = 1:numel(physioType)
    dtrig = diff(param.(physioType{p}).trig.data_filt);
    param.(physioType{p}).trig.average_rate = 1000 / mean(dtrig(abs(dtrig - mean(dtrig)) < 3 * std(dtrig))); 
end

fid = fopen(fullfile(datadir, 'physio.json'), 'w');
fprintf(fid, jsonencode(param, 'PrettyPrint', true));
fclose(fid);

%% END