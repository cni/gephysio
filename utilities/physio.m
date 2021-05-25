%% We downloaded some PPG data from
%
%  From the iris project, part of the LEANEW1 group
%  resting_2 physio data
%

%% Here are the scitran commands for finding and downloading the data

% Brian to complete
% st = scitran('cni');
%

%%
chdir(fullfile(gpRootPath,'data','gephysio1'));

% These are the pulsimetry values measured every 10 ms
foo = readmatrix('PPGData_cni_epi_0220202113_25_28_493');

ppg = foo(:,2);

% We convert from 10ms to 1 ms values
t = ((1:numel(ppg))-1)*10;

figure; 

% Choose a range of times to plot
start = 10000;
stop  = 30000;
r = logical(start < t) & logical(t < stop);

plot(t(r),ppg(r));
xlabel('Time (ms)');
ylabel('PPG value');

%% These are times that the PPG detects as an event

trig = readmatrix('PPGTrig_cni_epi_0220202113_25_28_493');

% Trigger in milliseconds
trig = trig*10;

lst = logical(start < trig) & logical (trig < stop);

% lst = (trig > r(1)) & (trig < r(end));

% Convert to milliseconds
% val = trig(lst)*10;

hold on;
mx = max(ppg(r)); mn = min(ppg(r));
% for ii=1:numel(lst)
%     line([trig(lst),trig(lst)],[mn,mx],'Color','k');
% end

set(gca,'ylim',[mn mx]);
plot(trig(lst),ones(size(trig(lst)))*mx*0.9,'ko');

%%  Respiration data

foo = readmatrix('RESPData_cni_epi_0220202113_25_28_493');

resp = foo(:,2);

% The samples were at 40 ms, not 10 ms
t = ((1:numel(resp))-1)*40;

figure; 

r = 9000:10000;

plot(t(r),resp(r));

trig_resp = readmatrix('RESPTrig_cni_epi_0220202113_25_28_493');

lst = (trig_resp > r(1)) & (trig_resp < r(end));

val = trig_resp(lst)*40;

hold on;

plot(val,ones(size(val))*3000,'ko');

%% Write out simple summary of the 

%% END