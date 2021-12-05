% example of standard linear chirps
clear all; close all;
Hz = 100 ; % sampling rate
L = 6 ; % time duration in second

t = [1/Hz:1/Hz:L]' ;
x1 = exp(2*pi*1i*(4*t.^2)); % first component
x2 = exp(2*pi*1i*(-pi*t.^2 + (24+6*pi).*t)); % second component

x = x1 + x2; % input signal

tt = -L:1/Hz:L;

% different windows
alpha1 = 1;
h1 = exp(-pi*alpha1*tt.^2'); % window g_0
Dh1 = dwindow(h1);
DDh1 = dwindow(Dh1);

alpha2 = 1;
h2 = tt'.^2.*exp(-pi*alpha2*tt.^2'); % window g_2
Dh2 = dwindow(h2);
DDh2 = dwindow(Dh2);

alpha3 = 5;
h3 = tt'.^2.*exp(-pi*alpha3*tt.^2'); % window g_3
Dh3 = dwindow(h3);
DDh3 = dwindow(Dh3);

alpha4 = 5;
h4 = tt'.^4.*exp(-pi*alpha4*tt.^2'); % window g_4
Dh4 = dwindow(h4);
DDh4 = dwindow(Dh4);

alpha5 = 5;
h5 = tt'.^6.*exp(-pi*alpha5*tt.^2'); % window g_5
Dh5 = dwindow(h5);
DDh5 = dwindow(Dh5);

% run SCT for different windows
[tfc1, tfrtic, tcrtic, tfrsq1, tfrsqtic] = sqSTCT(x, 0, 0.5, 2/length(x), 1, h1, Dh1, DDh1);
[tfc2, ~, ~, tfrsq2, ~] = sqSTCT(x, 0, 0.5, 2/length(x), 1, h2, Dh2, DDh2);
[tfc3, ~, ~, tfrsq3, ~] = sqSTCT(x, 0, 0.5, 2/length(x), 1, h3, Dh3, DDh3);
[tfc4, ~, ~, tfrsq4, ~] = sqSTCT(x, 0, 0.5, 2/length(x), 1, h4, Dh4, DDh4);
[tfc5, ~, ~, tfrsq5, ~] = sqSTCT(x, 0, 0.5, 2/length(x), 1, h5, Dh5, DDh5);

% show results between 1~5 second to avoid boundary effects
t_idx = 100:500;
t_show = t(t_idx);

return;
% Then run the block you want:
%% plot the upsampled signal and their IFs
scrsz = get(0,'ScreenSize');

t_up = [1/500:1/500:6]' ;
y1 = exp(2*pi*1i*(4*t_up.^2));
y2 = exp(2*pi*1i*(-pi*t_up.^2 + (24+6*pi).*t_up));

y = y1 + y2;

t_up_idx = 500:2500;
t_up_show = t_up(t_up_idx);

figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)*2/3])
plot(t_up_show,real(y1(t_up_idx)), 'color', 'k', 'linewidth', 1);
xlabel('time (s)')
set(gca,'fontsize',70)

figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)*2/3])
plot(t_up_show,real(y2(t_up_idx)), 'color', 'k', 'linewidth', 1);
xlabel('time (s)')
set(gca,'fontsize',70)

figure()
plot(t_up_show,real(y1(t_up_idx)+y2(t_up_idx)), 'color', 'k', 'linewidth', 1);
xlabel('time (s)')
set(gca,'fontsize',30)

if1 = 8*t;
if2 = -2*pi*t+(24+6*pi);

figure()
plot(t_show,if1(t_idx),'LineWidth',2);
hold on
plot(t_show,if2(t_idx),'LineWidth',2);
xlabel('time (s)');
ylabel('frequency (Hz)')
ylim([0 50])
set(gca,'fontsize',30)



%% plot the frequency-chirp rate slice of CT and SCT with g_0 and g_2 at 2s and 3s
figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfc1(:,:,200))), 1); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate'); 

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfrsq1(:,:,200))), 1); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate');

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfc1(:,:,300))), 1); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate'); 

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfrsq1(:,:,300))), 1); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate');

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfc2(:,:,300))), 1); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate'); 

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfrsq2(:,:,300))), 1); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate'); 


%% plot the contrast of CT and SCT with g_0 and g_2 at 3s (crossover)
figure()
A1 = abs(squeeze(tfc1(:,73,300)));
plot(tcrtic*Hz^2,A1, 'color', 'k');
xlim([-tcrtic(end)*Hz^2 tcrtic(end)*Hz^2])
xlabel('chirp rate')
set(gca,'fontsize',30)

figure()
A2 = abs(squeeze(tfc2(:,73,300)));
plot(tcrtic*Hz^2,A2, 'color', 'k');
xlim([-tcrtic(end)*Hz^2 tcrtic(end)*Hz^2])
xlabel('chirp rate')
set(gca,'fontsize',30)

figure()
A3 = abs(squeeze(tfrsq1(:,73,300)));
plot(tcrtic*Hz^2,A3, 'color', 'k');
xlim([-tcrtic(end)*Hz^2 tcrtic(end)*Hz^2])
xlabel('chirp rate')
set(gca,'fontsize',30)

figure()
A4 = abs(squeeze(tfrsq2(:,73,300)));
plot(tcrtic*Hz^2,A4, 'color', 'k');
xlim([-tcrtic(end)*Hz^2 tcrtic(end)*Hz^2])
xlabel('chirp rate')
set(gca,'fontsize',30)

%% spectral clustering to extract IFs and ICs
sigma = 60;
q = 0.9995;
mat = tfrsq2(:,:,t_idx);
idx = find(abs(mat)>quantile(abs(mat(:)),q));
[I1,I2,I3] = ind2sub(size(mat),idx);
% scatter3((I3+99)/Hz,tfrtic(I2)*Hz,tcrtic(I1)*Hz^2)
W = gaufunc([I1,I2,I3],[I1,I2,I3],sigma);
W = (W+W')/2;
% imagesc(W)
D = sum(W,2);
[evecs, evals] = eigs(W, diag(D),11);
label = kmeans(evecs(:,1:2),2);
idx1 = find(label==1);
idx2 = find(label==2);


line1 = zeros(length(t_show),2);
line2 = zeros(length(t_show),2);
line1_int = zeros(length(t_show),2);
line2_int = zeros(length(t_show),2);
df = (tfrtic(2)-tfrtic(1))*Hz;
dc = (tcrtic(2)-tcrtic(1))*Hz^2;

for i = t_idx-99
    t_idx1 = find(I3(idx1)==i);
    [max1, loc1] = max(abs(mat(idx(idx1(t_idx1)))));
    line1_int(i,1) = I2(idx1(t_idx1(loc1)));
    line1_int(i,2) = I1(idx1(t_idx1(loc1)));
    line1(i,1) = tfrtic(line1_int(i,1))*Hz;
    line1(i,2) = tcrtic(line1_int(i,2))*Hz^2;
    
    t_idx2 = find(I3(idx2)==i);
    [max2, loc2] = max(abs(mat(idx(idx2(t_idx2)))));
    line2_int(i,1) = I2(idx2(t_idx2(loc2)));
    line2_int(i,2) = I1(idx2(t_idx2(loc2)));
    line2(i,1) = tfrtic(line2_int(i,1))*Hz;
    line2(i,2) = tcrtic(line2_int(i,2))*Hz^2;
end

comb = [line1;line2];
rept = [t_show;t_show];
clr1 = repmat([0 0 1],numel(t_show),1);
clr2 = repmat([1 0 0],numel(t_show),1);
clr = [clr1;clr2];
figure
scatter3(rept,comb(:,1),comb(:,2),20,1-[0.85 0.85 0.85]);
ylim([0 50])
zlim([-tcrtic(end)*Hz^2 tcrtic(end)*Hz^2])
view(30,30)
xlabel('time (s)');
ylabel('frequency (Hz)');
zlabel('chirp rate');
colormap(1-gray)
set(gca,'fontsize',20)

% pad ones to match the original signal length
block1 = ones(99,2);
block2 = ones(100,2);
line1 = [block1; line1; block2];
line2 = [block1; line2; block2];
line1_int = [block1; line1_int; block2];
line2_int = [block1; line2_int; block2];


    %% 3d plot of chirplet transform with g_0
QN = 5 ;
D = tfc1(:,:,t_idx);
thresh = quantile(abs(D(:)),0.9999);
D(find(abs(D) < thresh * (10-QN+1)/10)) = thresh * (10-QN)/10 ;
 
for jj = 1: QN
    idx = find(abs(D) <= thresh * (10-jj+1)/10 & abs(D) > thresh * (10-jj)/10 );
    [I1,I2,I3] = ind2sub(size(D),idx);
    scatter3((I3+99)/Hz,tfrtic(I2)*Hz,tcrtic(I1)*Hz^2,20, [1 1 1]*(jj-1)/8, 'filled');
    hold on
end
 
view(30,30)
ylim([0 50])
zlim([-16.55 16.55])
xlabel('time (s)');
ylabel('frequency (Hz)');
zlabel('chirp rate');
colormap(1-gray)
colorbar
set(gca,'fontsize',20)

    %% 3d plot of SCT with g_0
QN = 5 ;
D = tfrsq1(:,:,t_idx);
thresh = quantile(abs(D(:)),0.9999);
D(find(abs(D) < thresh * (10-QN+1)/10)) = thresh * (10-QN)/10 ;
 
for jj = 1: QN
    idx = find(abs(D) <= thresh * (10-jj+1)/10 & abs(D) > thresh * (10-jj)/10 );
    [I1,I2,I3] = ind2sub(size(D),idx);
    scatter3((I3+99)/Hz,tfrtic(I2)*Hz,tcrtic(I1)*Hz^2,20, [1 1 1]*(jj-1)/8, 'filled');
    hold on
end
 
view(30,30)
ylim([0 50])
zlim([-16.55 16.55])
xlabel('time (s)');
ylabel('frequency (Hz)');
zlabel('chirp rate');
colormap(1-gray)
colorbar
set(gca,'fontsize',20)

    %% 3d plot of chirplet transform with g_2
QN = 5 ;
D = tfc2(:,:,t_idx);
thresh = quantile(abs(D(:)),0.9999);
D(find(abs(D) < thresh * (10-QN+1)/10)) = thresh * (10-QN)/10 ;
 
for jj = 1: QN
    idx = find(abs(D) <= thresh * (10-jj+1)/10 & abs(D) > thresh * (10-jj)/10 );
    [I1,I2,I3] = ind2sub(size(D),idx);
    scatter3((I3+99)/Hz,tfrtic(I2)*Hz,tcrtic(I1)*Hz^2,20, [1 1 1]*(jj-1)/8, 'filled');
    hold on
end
 
view(30,30)
ylim([0 50])
zlim([-16.55 16.55])
xlabel('time (s)');
ylabel('frequency (Hz)');
zlabel('chirp rate');
colormap(1-gray)
colorbar
set(gca,'fontsize',20)

    %% 3d plot of SCT with g_2
QN = 5 ;
D = tfrsq2(:,:,t_idx);
thresh = quantile(abs(D(:)),0.9999);
D(find(abs(D) < thresh * (10-QN+1)/10)) = thresh * (10-QN)/10 ;

for jj = 1: QN
    idx = find(abs(D) <= thresh * (10-jj+1)/10 & abs(D) > thresh * (10-jj)/10 );
    [I1,I2,I3] = ind2sub(size(D),idx);
    scatter3((I3+99)/Hz,tfrtic(I2)*Hz,tcrtic(I1)*Hz^2,20, [1 1 1]*(jj-1)/8, 'filled');
    hold on
end
 
view(30,30)
ylim([0 50])
zlim([-16.55 16.55])
xlabel('time (s)');
ylabel('frequency (Hz)');
zlabel('chirp rate');
colormap(1-gray)
colorbar
set(gca,'fontsize',30)


%%  plot the frequency-chirp rate slice and contrast of CT with g_3, g_4 and g_5 at 3s
figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfc3(:,:,300))), 1); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate'); 
set(gca,'fontsize',30)

figure()
A5 = abs(squeeze(tfc3(:,73,300)));
plot(tcrtic*Hz^2,A5, 'color', 'k');
xlim([-tcrtic(end)*Hz^2 tcrtic(end)*Hz^2])
xlabel('chirp rate')
set(gca,'fontsize',30)

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfc4(:,:,300))), 1); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate'); 
set(gca,'fontsize',30)

figure()
A6 = abs(squeeze(tfc4(:,73,300)));
plot(tcrtic*Hz^2,A6, 'color', 'k');
xlim([-tcrtic(end)*Hz^2 tcrtic(end)*Hz^2])
xlabel('chirp rate')
set(gca,'fontsize',30)

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfc5(:,:,300))), 1); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate'); 
set(gca,'fontsize',30)

figure()
A7 = abs(squeeze(tfc5(:,73,300)));
plot(tcrtic*Hz^2,A7, 'color', 'k');
xlim([-tcrtic(end)*Hz^2 tcrtic(end)*Hz^2])
xlabel('chirp rate')
set(gca,'fontsize',30)

%%  plot the frequency-chirp rate slice and contrast of SCT with g_3, g_4 and g_5 at 3s
figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfrsq3(:,:,300))), 1); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate'); 

figure()
A5 = abs(squeeze(tfrsq3(:,73,300)));
plot(tcrtic*Hz^2,A5, 'color', 'k');
xlim([-tcrtic(end)*Hz^2 tcrtic(end)*Hz^2])
xlabel('chirp rate')
set(gca,'fontsize',30)

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfrsq4(:,:,300))), 1); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate'); 

figure()
A6 = abs(squeeze(tfrsq4(:,73,300)));
plot(tcrtic*Hz^2,A6, 'color', 'k');
xlim([-tcrtic(end)*Hz^2 tcrtic(end)*Hz^2])
xlabel('chirp rate')
set(gca,'fontsize',30)

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfrsq5(:,:,300))), 1); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate'); 

figure()
A7 = abs(squeeze(tfrsq5(:,73,300)));
plot(tcrtic*Hz^2,A7, 'color', 'k');
xlim([-tcrtic(end)*Hz^2 tcrtic(end)*Hz^2])
xlabel('chirp rate')
set(gca,'fontsize',30)

%%  inverse SCT
ts = 300; % time in samples
[invpts] = SCTinverse(x, 0, 0.5, 2/length(x), 1, h1, Dh1, DDh1, ts, line2_int(ts,1), line2_int(ts,2), 3, 3);

figure
scatter(tfrtic(invpts(:,1))*Hz,tcrtic(invpts(:,2))*Hz^2,[],invpts(:,3).^1, 'filled');
xlabel('frequency (Hz)');ylabel('chirp-rate');
xlim([0 50])
colorbar;
set(gca,'fontsize',30)

%% reconstruction from 2nd-order SST
scrsz = get(0,'ScreenSize');
[x_recon1] = recon_sqSTCT(x, 0, 0.5, 2/length(x), 1, h1, Dh1, DDh1, ...
    line1_int(:,1), line1_int(:,2), 5, 5);
[x_recon2] = recon_sqSTCT(x, 0, 0.5, 2/length(x), 1, h1, Dh1, DDh1, ...
    line2_int(:,1), line2_int(:,2), 5, 5);
df = (tfrtic(2)-tfrtic(1))*Hz;
x_recon1 = x_recon1/Hz*df;
x_recon2 = x_recon2/Hz*df;
% legend fontsize 30
figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)*2/3])
plot(t_show,real(x1(t_idx)), 'color', 'k','linewidth',1);
hold on
plot(t_show,real(x_recon1(t_idx)));
ylim([-1.5 1.5])
legend('original','recon');
xlabel('time (s)')
set(gca,'fontsize',70)

figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)*2/3])
plot(t_show, real(x1(t_idx))-real(x_recon1(t_idx)), 'color', 'k');
ylim([-1.5 1.5])
xlabel('time (s)')
set(gca,'fontsize',70)

figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)*2/3])
plot(t_show,real(x2(t_idx)), 'color', 'k');
hold on
plot(t_show,real(x_recon2(t_idx)));
ylim([-1.5 1.5])
legend('original','recon');
xlabel('time (s)')
set(gca,'fontsize',70)

figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)*2/3])
plot(t_show, real(x2(t_idx))-real(x_recon2(t_idx)), 'color', 'k');
ylim([-1.5 1.5])
xlabel('time (s)')
set(gca,'fontsize',70)

S1_idx = 251:350;
errorS1 = norm(real(x1(S1_idx))-real(x_recon1(S1_idx)))/norm(real(x1(S1_idx)));
S2_idx = setdiff(t_idx,S1_idx);
errorS2 = norm(real(x1(S2_idx))-real(x_recon1(S2_idx)))/norm(real(x1(S2_idx)));

%% reconstruction by the group scheme from CT
scrsz = get(0,'ScreenSize');
chirp1 = line1(:,2);
chirp2 = line2(:,2);
if1 = line1(:,1);
if2 = line2(:,1);

rrcon = zeros(2,length(x));
scale = alpha1;
stctmat = tfc1;
for char = 1:length(x)
    tmp = zeros(2,2);
    tmp(1,1) = 1/sqrt(scale+1i*(tcrtic(line1_int(char,2))*Hz^2 - chirp1(char)))*exp(-pi*(tfrtic(line1_int(char,1))*Hz-if1(char))^2/(scale+1i*(tcrtic(line1_int(char,2))*Hz^2 - chirp1(char))));
    tmp(1,2) = 1/sqrt(scale+1i*(tcrtic(line1_int(char,2))*Hz^2 - chirp2(char)))*exp(-pi*(tfrtic(line1_int(char,1))*Hz-if2(char))^2/(scale+1i*(tcrtic(line1_int(char,2))*Hz^2 - chirp2(char))));
    tmp(2,1) = 1/sqrt(scale+1i*(tcrtic(line2_int(char,2))*Hz^2 - chirp1(char)))*exp(-pi*(tfrtic(line2_int(char,1))*Hz-if1(char))^2/(scale+1i*(tcrtic(line2_int(char,2))*Hz^2 - chirp1(char))));
    tmp(2,2) = 1/sqrt(scale+1i*(tcrtic(line2_int(char,2))*Hz^2 - chirp2(char)))*exp(-pi*(tfrtic(line2_int(char,1))*Hz-if2(char))^2/(scale+1i*(tcrtic(line2_int(char,2))*Hz^2 - chirp2(char))));
    
    xtmp = [stctmat(line1_int(char,2),line1_int(char,1),char); stctmat(line2_int(char,2),line2_int(char,1),char)]/Hz;
    rrcon(:,char) = tmp\xtmp;
end

figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)*2/3])
plot(t_show,real(x1(t_idx)), 'color', 'k');
xlabel('time (s)');
hold on
plot(t_show,real(rrcon(1,t_idx)));
ylim([-1.5 1.5])
xlabel('time (s)');
legend('original','recon');
set(gca,'fontsize',70)

figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)*2/3])
plot(t_show, real(x1(t_idx)')-real(rrcon(1,t_idx)), 'color', 'k');
ylim([-1.5 1.5])
xlabel('time (s)');
set(gca,'fontsize',70)

figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)*2/3])
plot(t_show,real(x2(t_idx)), 'color', 'k');
xlabel('time (s)');
hold on
plot(t_show,real(rrcon(2,t_idx)));
ylim([-1.5 1.5])
legend('original','recon');
xlabel('time (s)')
set(gca,'fontsize',70)

figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)*2/3])
plot(t_show, real(x2(t_idx)')-real(rrcon(2,t_idx)), 'color', 'k');
ylim([-1.5 1.5])
xlabel('time (s)');
set(gca,'fontsize',70)

G1_idx = 251:350;
errorG1 = norm(real(x1(G1_idx)')-real(rrcon(1,G1_idx)))/norm(real(x1(G1_idx)));
G2_idx = setdiff(t_idx,G1_idx);
errorG2 = norm(real(x1(G2_idx)')-real(rrcon(1,G2_idx)))/norm(real(x1(G2_idx)));

    %% reassignment
[tfr,rtfr,hat] = tfrrsp(x,1:length(x),length(x)/2,h1,0);
rtfr = rtfr(1:length(x)/4,:);

imageSQ(t(t_idx), 0:Hz*2/length(x):(Hz/2-Hz*2/length(x)), abs(rtfr(:,t_idx)),0.9999);axis xy; colormap(1-gray) ;
xlabel('time (s)'); ylabel('frequency (Hz)');

    %% 1st-order SST
[tfr0, tfrtic0, tfrsq0, tfrsqtic0] = sqSTFT(x, 0, 0.5, 1/length(x), 1, h1, Dh1);
figure()
imageSQ(t(t_idx), Hz*tfrsqtic0, abs(tfrsq0(:,t_idx)), 0.9999); axis xy; colormap(1-gray); 
xlabel('time (s)'); ylabel('frequency (Hz)');

    %% 2nd-order SST
[tfr0, tfrtic0, tfrsq0, tfrsq2nd0, tfrsqtic0] = sqSTFTbase2nd(x, 0, 0.5, 2/length(x), 1, h1, Dh1, DDh1);
figure()
imageSQ(t(t_idx), Hz*tfrsqtic0, abs(tfrsq2nd0(:,t_idx)), 0.9999); axis xy; colormap(1-gray); 
xlabel('time (s)'); ylabel('frequency (Hz)'); 

    %% projection of CT onto TF plane
tfproj = zeros(size(tfc1,2),size(tfc1,3));
for i = 1:size(tfproj,1)
    for j = 1:size(tfproj,2)
        tfproj(i,j) = sum(abs(squeeze(tfc1(:,i,j))));
    end
end
figure()
imageSQ(t(t_idx), Hz*tfrtic, abs(tfproj(:,t_idx)), 0.9999); axis xy; colormap(1-gray); 
xlabel('time (s)'); ylabel('frequency (Hz)'); 

    %% projection of SCT onto TF plane
tfproj = zeros(size(tfrsq1,2),size(tfrsq1,3));
for i = 1:size(tfproj,1)
    for j = 1:size(tfproj,2)
        tfproj(i,j) = sum(abs(squeeze(tfrsq1(:,i,j))));
    end
end
figure()
imageSQ(t(t_idx), Hz*tfrtic, abs(tfproj(:,t_idx)), 0.9999); axis xy; colormap(1-gray); 
xlabel('time (s)'); ylabel('frequency (Hz)'); 

%% bipartite extraction (not used)
% line1 = zeros(length(t),2);
% line2 = zeros(length(t),2);
% line1_int = zeros(length(t),2);
% line2_int = zeros(length(t),2);
% C = tfrsq2; % for curve extraction
% %C(abs(C)<=quantile(abs(mat(:)),q))=0;
% cth = ceil(size(C,1)/2);
% C1 = C;
% C2 = C;
% C1(1:cth,:,:) = 0;
% C2(cth+1:end,:,:) = 0;
% df = (tfrtic(2)-tfrtic(1))*Hz;
% dc = (tcrtic(2)-tcrtic(1))*Hz^2;
% 
% for i = 1:length(t)
%     E1 = abs(squeeze(C1(:,:,i)));
%     [max1, idx1] = max(E1(:));
%     freq1 = floor(idx1/size(C,1))+1;
%     chirp1 = mod(idx1,size(C,1));
%     line1_int(i,1) = freq1;
%     line1_int(i,2) = chirp1;
%     freq1 = tfrtic(freq1);
%     chirp1 = tcrtic(chirp1);
%     line1(i,1) = freq1*Hz;
%     line1(i,2) = chirp1*Hz^2;
%     
%     E2 = abs(squeeze(C2(:,:,i)));
%     [max2, idx2] = max(E2(:));
%     freq2 = floor(idx2/size(C,1))+1;
%     chirp2 = mod(idx2,size(C,1));
%     line2_int(i,1) = freq2;
%     line2_int(i,2) = chirp2;
%     freq2 = tfrtic(freq2)*Hz;
%     chirp2 = tcrtic(chirp2)*Hz^2;
%     line2(i,1) = freq2;
%     line2(i,2) = chirp2;
% end
% 
% comb = [line1(t_idx,:);line2(t_idx,:)];
% rept = [t_show;t_show];
% clr1 = repmat([0 0 1],numel(t_show),1);
% clr2 = repmat([1 0 0],numel(t_show),1);
% clr = [clr1;clr2];
% figure
% scatter3(rept,comb(:,1),comb(:,2),20,clr);
% ylim([0 50])
% zlim([-tcrtic(end)*Hz^2 tcrtic(end)*Hz^2])
% view(30,30)
% xlabel('time (s)');
% ylabel('frequency (Hz)');
% zlabel('chirp rate');
% set(gca,'fontsize',20)

%% generate animation of rotation
filename = 'rotation0.gif';

QN = 5 ;
D = tfc1(:,:,t_idx);
thresh = quantile(abs(D(:)),0.9999);
D(find(abs(D) < thresh * (10-QN+1)/10)) = thresh * (10-QN)/10 ;
 
idx = find(abs(D) <= thresh * (10-1+1)/10 & abs(D) > thresh * (10-1)/10 );
[I11,I12,I13] = ind2sub(size(D),idx);

idx = find(abs(D) <= thresh * (10-2+1)/10 & abs(D) > thresh * (10-2)/10 );
[I21,I22,I23] = ind2sub(size(D),idx);

idx = find(abs(D) <= thresh * (10-3+1)/10 & abs(D) > thresh * (10-3)/10 );
[I31,I32,I33] = ind2sub(size(D),idx);

idx = find(abs(D) <= thresh * (10-4+1)/10 & abs(D) > thresh * (10-4)/10 );
[I41,I42,I43] = ind2sub(size(D),idx);

idx = find(abs(D) <= thresh * (10-5+1)/10 & abs(D) > thresh * (10-5)/10 );
[I51,I52,I53] = ind2sub(size(D),idx);

h = figure;
angle = 0:3.6:360;
for i = 1:length(angle)
    clf;
    scatter3((I13+99)/Hz,tfrtic(I12)*Hz,tcrtic(I11)*Hz^2,20, [1 1 1]*(1-1)/8, 'filled');
    hold on
    scatter3((I23+99)/Hz,tfrtic(I22)*Hz,tcrtic(I21)*Hz^2,20, [1 1 1]*(2-1)/8, 'filled');
    scatter3((I33+99)/Hz,tfrtic(I32)*Hz,tcrtic(I31)*Hz^2,20, [1 1 1]*(3-1)/8, 'filled');
    scatter3((I43+99)/Hz,tfrtic(I42)*Hz,tcrtic(I41)*Hz^2,20, [1 1 1]*(4-1)/8, 'filled');
    scatter3((I53+99)/Hz,tfrtic(I52)*Hz,tcrtic(I51)*Hz^2,20, [1 1 1]*(5-1)/8, 'filled');
    view(angle(i),30)
    ylim([0 50])
    zlim([-16.55 16.55])
    xlabel('time (s)');
    ylabel('frequency (Hz)');
    zlabel('chirp rate');
    set(gca,'fontsize',15)
    drawnow
    
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    % Write to the GIF File
    if i == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end
