foo = readmatrix('PPGData_cni_epi_0220202113_25_28_493');

ppg = foo(:,2);

t = ((1:numel(ppg))-1)*10;

% ieNewGraphWin;

figure; 

r = 9000:10000;

plot(t(r),ppg(r));

trig = readmatrix('PPGTrig_cni_epi_0220202113_25_28_493');

lst = (trig > r(1)) & (trig < r(end));

val = trig(lst);

hold on;

plot(val*10,ones(size(val))*100,'ko');

%%
foo = readmatrix('RESPData_cni_epi_0220202113_25_28_493');

resp = foo(:,2);

t = ((1:numel(resp))-1)*40;

% ieNewGraphWin;

figure; 

r = 9000:10000;

plot(t(r),resp(r));

trig_resp = readmatrix('RESPTrig_cni_epi_0220202113_25_28_493');

lst = (trig_resp > r(1)) & (trig_resp < r(end));

val = trig_resp(lst);

hold on;

plot(val*10,ones(size(val))*3000,'ko');

