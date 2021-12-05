% example of a signal with time-varying amplitudes and chirp rates
clear all; close all ;  

	% the sampling rate for the simulated signal
Hz = 100 ;
	% the sampling time of the simulated signal
time = [1/Hz:1/Hz:10]' ;
	% the number of sampling points of the simulated signal
N = length(time) ;

    % fix the random seed for the reproducibility issue
initstate(1) ;


	% the amplitude modulation of the simulated signal
	% simulate 2 oscillatory components with dynamics
am1 = smooth(cumsum(randn(N,1)) ./ Hz, 200, 'loess') ;
am1 = 2 + am1 ./ max(abs(am1)) ;
am2 = smooth(cumsum(randn(N,1)) ./ Hz, 200, 'loess') ;
am2 = 2 + am2 ./ max(abs(am2)) ;


    %% the instantaneous chirp-rate of the simulated signal
chirp1 = smooth(cumsum(randn(N,1)) ./ Hz, 400, 'loess') ;
chirp1 = 4.5 + chirp1 ./ max(abs(chirp1))/5 ;
chirp2 = smooth(cumsum(randn(N,1)) ./ Hz, 300, 'loess') ;
chirp2 = -4 + chirp2 ./ max(abs(chirp2))/4 ;
 
    %% the instantaneous frequency of the simulated signal
if1 = cumsum(chirp1) / Hz + 2; 
if2 = cumsum(chirp2) / Hz + 44; 
phi1 = cumsum(if1) / Hz ; 
phi2 = cumsum(if2) / Hz ; 
 
	%% the simulated signal.
s1 = am1 .* exp(2*pi*1i*phi1) ; 
s2 = am2 .* exp(2*pi*1i*phi2) ; 
clean = s1 + s2 ;

%%
	% add noise
sigma = 1 ;%sqrt( var(clean)*10.^( -snrdb /10 ) );
noise = random('T',4,N,1) ;
noise = sigma * noise ; 
var(noise)
snrdb = 20 * log10(std(clean)./std(noise)) ;
fprintf(['snrdb = ',num2str(snrdb),'\n']) ;

	% simulated observed time series
xm = clean + noise ; % input signal

% show results between 1~9 second to avoid boundary effects
t_idx = 100:900;
t_show = time(t_idx);

    %% plot signals and their IFs and ICs
scrsz = get(0,'ScreenSize');

figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)])
plot(t_show,real(s1(t_idx)), 'color', 'k');
xlabel('time (s)');
xlim([1 9]);
set(gca,'fontsize',100)

figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)])
plot(t_show,real(s2(t_idx)), 'color', 'k');
xlabel('time (s)');
xlim([1 9]);
set(gca,'fontsize',100)

figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)])
plot(t_show,real(s1(t_idx)+s2(t_idx)), 'color', 'k');
xlabel('time (s)');
xlim([1 9]);
ylim([-5 5]);
set(gca,'fontsize',100)

figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)])
plot(t_show,real(xm(t_idx)), 'color', 'k');
xlabel('time (s)');
xlim([1 9]);
ylim([-5 5]);
set(gca,'fontsize',100)

figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)])
plot(t_show,am1(t_idx), 'color', 'b', 'linewidth', 1);
hold on;
plot(t_show,am2(t_idx), 'color', 'r', 'linewidth', 1);
xlim([1 9]);
set(gca,'fontsize',100)
    
figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)])
plot(t_show,if1(t_idx), 'color', 'b', 'linewidth', 1);
ylabel('frequency (Hz)');
hold on;
plot(t_show,if2(t_idx), 'color', 'r', 'linewidth', 1);
xlim([1 9]);
ylim([0 50]);
xlabel('time (s)');
ylabel('frequency (Hz)');
set(gca,'fontsize',100)

figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)*2/3])
plot(t_show,chirp1(t_idx), 'color', 'b', 'linewidth', 1);
hold on;
plot(t_show,chirp2(t_idx), 'color', 'r', 'linewidth', 1);
xlim([1 9]);
xlabel('time (s)');
ylabel('chirp rate');
xlim([1 9])
set(gca,'fontsize',100)

    %% run the SCT with g_0 and g_2
alpha1 = 1;
tt = -8:1/Hz:8;
h1 = exp(-pi*alpha1*tt.^2'); % window g_0
Dh1 = dwindow(h1); 
DDh1 = dwindow(Dh1);

alpha2 = 1;
tt = -8:1/Hz:8;
h2 = tt'.^2.*exp(-pi*alpha2*tt.^2'); % window g_2
Dh2 = dwindow(h2); 
DDh2 = dwindow(Dh2);

[tfc1, tfrtic, tcrtic, tfrsq1, ~] = sqSTCT(xm, 0, 0.5, 2.5/length(xm), 1, h1, Dh1, DDh1);
[tfc2, ~, ~, tfrsq2, ~] = sqSTCT(xm, 0, 0.5, 2.5/length(xm), 1, h2, Dh2, DDh2);
return;
%% frequency-chirp rate slice of CT and SCT with g_0 and g_2 at crossover time and the contrast
figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfc1(:,:,502))), 1); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate'); 

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfrsq1(:,:,502))), 1); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate'); 

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfc2(:,:,502))), 1); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate');

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfrsq2(:,:,502))), 1); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate');

figure()
plot(tcrtic*Hz^2,abs(squeeze(tfrsq1(:,98,502))), 'color', 'k')
xlabel('chirp rate')
set(gca,'fontsize',30)

figure()
plot(tcrtic*Hz^2,abs(squeeze(tfrsq2(:,98,502))), 'color', 'k')
xlabel('chirp rate')
set(gca,'fontsize',30)

    %% time-frequency-chirp_rate curve (ground truth, unit in sample point)
df = (tfrtic(2)-tfrtic(1))*Hz;
dc = (tcrtic(2)-tcrtic(1))*Hz^2;
c0 = tcrtic(1)*Hz^2;
line1 = zeros(length(time),2);
line2 = zeros(length(time),2);

for j = 1:length(time)
    line1(j,1) = round(if1(j)/df)+1;
    line1(j,2) = round((chirp1(j)-c0)/dc)+1;
    
    line2(j,1) = round(if2(j)/df)+1;
    line2(j,2) = round((chirp2(j)-c0)/dc)+1;
end

    %% TFC ridges from the SCT with g_2
t = time;
line1 = zeros(length(t),2);
line2 = zeros(length(t),2);
line1_int = zeros(length(t),2);
line2_int = zeros(length(t),2);
C = tfrsq2; % for curve extraction
cth = ceil(size(C,1)/2);
C1 = C;
C2 = C;
C1(1:cth,:,:) = 0;
C2(cth+1:end,:,:) = 0;
df = (tfrtic(2)-tfrtic(1))*Hz;
dc = (tcrtic(2)-tcrtic(1))*Hz^2;

for i = 1:length(t)
    E1 = abs(squeeze(C1(:,:,i)));
    [max1, idx1] = max(E1(:));
    freq1 = floor(idx1/size(C,1))+1;
    chirp1 = mod(idx1,size(C,1));
    line1_int(i,1) = freq1;
    line1_int(i,2) = chirp1;
    freq1 = tfrtic(freq1);
    chirp1 = tcrtic(chirp1);
    line1(i,1) = freq1*Hz;
    line1(i,2) = chirp1*Hz^2;
    
    E2 = abs(squeeze(C2(:,:,i)));
    [max2, idx2] = max(E2(:));
    freq2 = floor(idx2/size(C,1))+1;
    chirp2 = mod(idx2,size(C,1));
    line2_int(i,1) = freq2;
    line2_int(i,2) = chirp2;
    freq2 = tfrtic(freq2)*Hz;
    chirp2 = tcrtic(chirp2)*Hz^2;
    line2(i,1) = freq2;
    line2(i,2) = chirp2;
end

comb = [line1(t_idx,:);line2(t_idx,:)];
rept = [t_show;t_show];
figure
scatter3(rept,comb(:,1),comb(:,2),20,comb(:,2));
view(30,30)
xlabel('time (s)');
ylabel('frequency (Hz)');
zlabel('chirp rate');
set(gca,'fontsize',20)

chirp1 = line1(:,2);
chirp2 = line2(:,2);
if1 = line1(:,1);
if2 = line2(:,1);
line1 = line1_int;
line2 = line2_int;

    %% reconstruction (run one of the previous two blocks before this block)
scrsz = get(0,'ScreenSize');
rrcon = zeros(2,length(clean));
scale = alpha1;
stctmat = tfc1;
for char = 1:length(clean)
    tmp = zeros(2,2);
    tmp(1,1) = 1/sqrt(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirp1(char)))*exp(-pi*(tfrtic(line1(char,1))*Hz-if1(char))^2/(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirp1(char))));
    tmp(1,2) = 1/sqrt(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirp2(char)))*exp(-pi*(tfrtic(line1(char,1))*Hz-if2(char))^2/(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirp2(char))));
    tmp(2,1) = 1/sqrt(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirp1(char)))*exp(-pi*(tfrtic(line2(char,1))*Hz-if1(char))^2/(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirp1(char))));
    tmp(2,2) = 1/sqrt(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirp2(char)))*exp(-pi*(tfrtic(line2(char,1))*Hz-if2(char))^2/(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirp2(char))));
    
    xtmp = [stctmat(line1(char,2),line1(char,1),char); stctmat(line2(char,2),line2(char,1),char)]/Hz;
    rrcon(:,char) = tmp\xtmp;
end

figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)*2/3])
plot(t_show,real(s1(t_idx)), 'color', 'k');
xlabel('time (s)');
hold on
plot(t_show,real(rrcon(1,t_idx)));
xlabel('time (s)');
xlim([1 9])
ylim([-3 3])
legend('original','recon');
set(gca,'fontsize',70)

figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)*2/3])
plot(t_show, real(s1(t_idx)')-real(rrcon(1,t_idx)), 'color', 'k');
xlim([1 9])
ylim([-3 3])
xlabel('time (s)');
set(gca,'fontsize',70)

figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)*2/3])
plot(t_show,real(s2(t_idx)), 'color', 'k');
xlabel('time (s)');
hold on
plot(t_show,real(rrcon(2,t_idx)));
legend('original','recon');
xlabel('time (s)')
xlim([1 9])
ylim([-3 3])
set(gca,'fontsize',70)

figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)*2/3])
plot(t_show, real(s2(t_idx)')-real(rrcon(2,t_idx)), 'color', 'k');
xlim([1 9])
ylim([-3 3])
xlabel('time (s)');
set(gca,'fontsize',70)

    %% 2nd-order SST
[tfr0, tfrtic0, tfrsq0, tfrsq2nd0, tfrsqtic0] = sqSTFTbase2nd(xm, 0, 0.5, 2.5/length(xm), 1, h1, Dh1, DDh1);
figure()
imageSQ(time(t_idx), Hz*tfrsqtic0, abs(tfrsq2nd0(:,t_idx)), 0.9999); axis xy; colormap(1-gray); 
xlabel('time (s)'); ylabel('frequency (Hz)'); 

    %% projection of SCT onto TF plane
tfproj = zeros(size(tfrsq1,2),size(tfrsq1,3));
for i = 1:size(tfproj,1)
    for j = 1:size(tfproj,2)
        tfproj(i,j) = sum(abs(squeeze(tfrsq1(:,i,j))));
    end
end
figure()
imageSQ(time(t_idx), Hz*tfrtic, abs(tfproj(:,t_idx)), 0.9999); axis xy; colormap(1-gray); 
xlabel('time (s)'); ylabel('frequency (Hz)'); 
