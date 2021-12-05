% example of a wolf howling signal
clear all; close all;

[x,Hz] = audioread('Wolf.wav');
x = downsample(x,8); % downsample the signal by 8
Hz = Hz/8;
samples = [round(16*Hz):round(17*Hz)-1];
xm = x(samples);

hop = 1;
L = 8;
tt = -L:1/20:L;

alpha1 = 0.1; 
h1 = exp(-pi*alpha1*tt.^2'); % window g_0
Dh1 = dwindow(h1);
DDh1 = dwindow(Dh1);

alpha2 = 0.1;
h2 = tt'.^2.*exp(-pi*alpha2*tt.^2'); % window g_2
Dh2 = dwindow(h2);
DDh2 = dwindow(Dh2);


    %% run the 2nd SST, SCT
[tfr0, tfrtic0, tfrsq0, tfrsq2nd, tfrsqtic0] = sqSTFTbase2nd(xm, 0, 0.5, 2/length(xm), hop, h1, Dh1, DDh1);
[tfc1, tfrtic, tcrtic, tfrsq1, tfrsqtic] = sqSTCT(xm, 0, 0.5, 2/length(xm), hop, h1, Dh1, DDh1);
[tfc2, ~, ~, tfrsq2, ~] = sqSTCT(xm, 0, 0.5, 2/length(xm), hop, h2, Dh2, DDh2);
df = (tfrtic(2)-tfrtic(1))*Hz;
dc = (tcrtic(2)-tcrtic(1))*Hz^2;
time_stamp = hop/Hz;
return;

% Then run the block you want:
    %% figure 12
figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfc1(:,:,161))).^2, 0.9999); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate');

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfrsq1(:,:,161))).^2, 0.9999); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate');

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfc2(:,:,161))).^2, 0.9999); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate');

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfrsq2(:,:,161))).^2, 0.9999); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate');

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfc1(:,:,431))).^2, 0.9999); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate');

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfrsq1(:,:,431))).^2, 0.9999); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate');

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfc2(:,:,431))).^2, 0.9999); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate');

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfrsq2(:,:,431))).^2, 0.9999); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate');

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfc2(:,:,571))).^2, 0.9999); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate');

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfrsq2(:,:,571))).^2, 0.9999); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate');


    %% plot the STFT and the 2nd-order SST
figure() % STFT
imageSQ(16:time_stamp:16+time_stamp*(size(tfrsq0,2)-1), Hz*tfrtic0, abs(tfr0), 0.9999); axis xy; colormap(1-gray); 
xlabel('time (s)');ylabel('frequency (Hz)');

figure() % 2nd-order SST
imageSQ(16:time_stamp:16+time_stamp*(size(tfrsq0,2)-1), Hz*tfrtic0, abs(tfrsq2nd), 0.9999); axis xy; colormap(1-gray); 
xlabel('time (s)');ylabel('frequency (Hz)');

    %% plot the signal
scrsz = get(0,'ScreenSize');

figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(4)/3])
plot(16:1/Hz:16+(length(xm)-1)/Hz,xm, 'color', 'k', 'linewidth', 1);
xlabel('time (s)')
set(gca,'fontsize',30)



    %% 3d plot of chirplet transform with g_0
QN = 6 ;
D = tfc2;
thresh = quantile(abs(D(:)),0.9999);
D(find(abs(D) < thresh * (10-QN+1)/10)) = thresh * (10-QN)/10 ;
 
for jj = 1: QN
    idx = find(abs(D) <= thresh * (10-jj+1)/10 & abs(D) > thresh * (10-jj)/10 );
    [I1,I2,I3] = ind2sub(size(D),idx);
    scatter3(16+(I3-1)/Hz,tfrtic(I2)*Hz,tcrtic(I1)*Hz^2,20, [1 1 1]*(jj-1)/8, 'filled');
    hold on
end
 
view(30,30)
ylim([0 500])
xlabel('time (s)');
ylabel('frequency (Hz)');
zlabel('chirp rate');
colormap(1-gray)
colorbar
set(gca,'fontsize',20)

    %% projection of CT onto TF plane
tfproj = zeros(size(tfc1,2),size(tfc1,3));
for i = 1:size(tfproj,1)
    for j = 1:size(tfproj,2)
        tfproj(i,j) = sum(abs(squeeze(tfc1(:,i,j))));
    end
end
figure()
imageSQ(16:1/Hz:17-1/Hz, Hz*tfrtic, abs(tfproj), 0.9999); axis xy; colormap(1-gray); 
xlabel('time (s)'); ylabel('frequency (Hz)'); 

    %% projection of SCT onto TF plane
tfproj = zeros(size(tfrsq1,2),size(tfrsq1,3));
for i = 1:size(tfproj,1)
    for j = 1:size(tfproj,2)
        tfproj(i,j) = sum(abs(squeeze(tfrsq1(:,i,j))));
    end
end
figure()
imageSQ(16:1/Hz:17-1/Hz, Hz*tfrtic, abs(tfproj), 1); axis xy; colormap(1-gray); 
xlabel('time (s)'); ylabel('frequency (Hz)'); 

    %% scatter plot of SCT
val = 0.9995;
thresh = quantile(abs(tfrsq1(:)),val);
idx = find(tfrsq1>thresh);
[I1,I2,I3] = ind2sub(size(tfrsq1),idx);
scatter3(16+(I3-1)/Hz,tfrtic(I2)*Hz,tcrtic(I1)*Hz^2,20,abs(tfrsq1(idx)));
view(30,30)
ylim([0 500])
xlabel('time (s)');
ylabel('frequency (Hz)');
zlabel('chirp rate');
set(gca,'fontsize',20)
