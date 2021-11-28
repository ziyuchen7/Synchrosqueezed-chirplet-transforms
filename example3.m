% example of signals with parabolic IFs
clear all; close all;
Hz = 100 ;
L = 6 ;

t = [1/Hz:1/Hz:L]' ;
phi1 = 8*t + (t-3).^3/3;
phi2 = 10*t - (t-3).^3/3;
x1 = exp(2*pi*1i*phi1);
x2 = exp(2*pi*1i*phi2);

x = x1 + x2;

tt = -L:1/Hz:L;

alpha1 = 5;
h1 = exp(-pi*alpha1*tt.^2'); % window g_0
Dh1 = dwindow(h1);
DDh1 = dwindow(Dh1);

alpha2 = 5;
h2 = tt'.^2.*exp(-pi*alpha2*tt.^2'); % window g_2
Dh2 = dwindow(h2);
DDh2 = dwindow(Dh2);


[tfc1, tfrtic, tcrtic, tfrsq1, tfrsqtic] = sqSTCT(x, 0, 0.5, 2/length(x), 1, h1, Dh1, DDh1);
[tfc2, ~, ~, tfrsq2, ~] = sqSTCT(x, 0, 0.5, 2/length(x), 1, h2, Dh2, DDh2);
return;
t_idx = 100:500;
t_show = t(t_idx);

%% figure 12

if1 = 8 + (t-3).^2;
if2 = 10 - (t-3).^2;

figure()
plot(t_show,if1(t_idx));
hold on
plot(t_show,if2(t_idx));
xlabel('time (s)');
ylabel('frequency (Hz)')
set(gca,'fontsize',30)

chirp1 = 2*(t-3);
chirp2 = -2*(t-3);

figure()
plot(t_show,chirp1(t_idx));
hold on
plot(t_show,chirp2(t_idx));
xlabel('time (s)');
ylabel('chirp rate')
set(gca,'fontsize',30)

%% figure 12
figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfrsq1(:,:,100))), 1); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate');

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfrsq1(:,:,200))), 1); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate');

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfrsq1(:,:,300))), 1); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate');

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfrsq2(:,:,100))), 1); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate');

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfrsq2(:,:,200))), 1); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate');

figure()
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfrsq2(:,:,300))), 1); axis xy; colormap(1-gray); 
xlabel('frequency (Hz)'); ylabel('chirp rate');


%% 
figure()
A1 = abs(squeeze(tfc1(:,73,300)));
plot(tcrtic*Hz^2,A1);
xlabel('chirp rate')
set(gca,'fontsize',20)

figure()
A2 = abs(squeeze(tfc2(:,73,300)));
plot(tcrtic*Hz^2,A2);
xlabel('chirp rate')
set(gca,'fontsize',20)

figure()
A3 = abs(squeeze(tfrsq1(:,73,300)));
plot(tcrtic*Hz^2,A3);
xlabel('chirp rate')
set(gca,'fontsize',20)

figure()
A4 = abs(squeeze(tfrsq2(:,73,300)));
plot(tcrtic*Hz^2,A4);
xlabel('chirp rate')
set(gca,'fontsize',20)

