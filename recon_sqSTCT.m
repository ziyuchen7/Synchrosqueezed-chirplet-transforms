function [x_recon] = recon_sqSTCT(x, lowFreq, highFreq, alpha, tDS, h, Dh, DDh, freq_idx, chirp_idx, fband, cband)

x_recon = zeros(size(x));
[xrow,xcol] = size(x) ;
t = [1:length(x)] ;
tLen = length(t(1:tDS:length(x))) ;

	% for tfr
N = length([-0.5+alpha:alpha:0.5]) ;
crate = ([1:N-1]-ceil(N/2))/N^2;

	% for tfrsq
Lidx = ceil( (N/2)*(lowFreq/0.5) ) + 1 ; 
Hidx = floor( (N/2)*(highFreq/0.5) ) ;  
fLen = Hidx - Lidx + 1 ;
%====================================================================
	%% check input signals
if (xcol~=1)
    error('X must have only one column');
elseif highFreq > 0.5
    error('TopFreq must be a value in [0, 0.5]');
elseif (tDS < 1) || (rem(tDS,1)) 
    error('tDS must be an integer value >= 1');
end

[hrow,hcol] = size(h); Lh = (hrow-1)/2; 
if (hcol~=1)||(rem(hrow,2)==0)
    error('H must be a smoothing window with odd length');
end
ht = [-Lh:Lh] ;

%====================================================================
	%% run STFT and reassignment rule
Ex = mean(abs(x(min(t):max(t))).^2);
Threshold = 1.0e-8*Ex;  % originally it was 1e-6*Ex

for tidx = 1:tLen
    chirp = crate(chirp_idx(tidx));
    
    % ti is the current time
    ti = t((tidx-1)*tDS+1);
    % tau is the relevant index associated with ti
    tau = -min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);
    % indices is the absolute index in the "evaluation window"
    indices= rem(N+tau,N)+1;
    norm_h = norm(h(Lh+1+tau));
    
    tf0 = zeros(N, 1) ; tf1 = zeros(N, 1) ; tf2 = zeros(N, 1) ;
    tfx0 = zeros(N, 1) ; tfx1 = zeros(N, 1) ; tfx2 = zeros(N,1);
    
    tf0(indices) = x(ti+tau).*conj(h(Lh+1+tau)).*exp(-pi*1i*chirp.*...
        (ht(Lh+1+tau)').^2);
    tf1(indices) = x(ti+tau).*conj(Dh(Lh+1+tau)).*exp(-pi*1i*chirp.*...
        (ht(Lh+1+tau)').^2);
    tf2(indices) = x(ti+tau).*conj(DDh(Lh+1+tau)).*exp(-pi*1i*chirp.*...
        (ht(Lh+1+tau)').^2);
    tfx0(indices) = x(ti+tau).*conj(h(Lh+1+tau)).*ht(Lh+1+tau)'.*...
        exp(-pi*1i*chirp.*(ht(Lh+1+tau)').^2);
    tfx1(indices) = x(ti+tau).*conj(Dh(Lh+1+tau)).*ht(Lh+1+tau)'.*...
        exp(-pi*1i*chirp.*(ht(Lh+1+tau)').^2);
    tfx2(indices) = x(ti+tau).*conj(h(Lh+1+tau)).*(ht(Lh+1+tau).^2)'.*...
        exp(-pi*1i*chirp.*(ht(Lh+1+tau)').^2);
    
    tf0 = fft(tf0) ; tf0 = tf0(1:N/2) ;
    tf1 = fft(tf1) ; tf1 = tf1(1:N/2) ;
    tf2 = fft(tf2) ; tf2 = tf2(1:N/2) ;
    tfx0 = fft(tfx0) ; tfx0 = tfx0(1:N/2) ;
    tfx1 = fft(tfx1) ; tfx1 = tfx1(1:N/2) ;
    tfx2 = fft(tfx2) ; tfx2 = tfx2(1:N/2) ;
    
    % get the first order omega
    lambda0 = (tf0.*tf2 - 4*pi*1i*chirp*tf0.*tfx1 - 2*pi*1i*chirp*...
        tf0.*tf0 + (2*pi*1i*chirp)^2*tf0.*tfx2 - tf1.*tf1 -(2*pi*1i*chirp)^2*...
        tfx0.*tfx0 + 4*pi*1i*chirp*tf1.*tfx0)./(-tf0.*tfx1 + 2*pi*1i*chirp*tf0.*tfx2 +...
        tfx0.*tf1 - 2*pi*1i*chirp*tfx0.*tfx0)./(2.0*pi);
    
    lambda = round(N^2 * imag(lambda0));
    
    omega = round(N * imag(tf1./tf0./(2.0*pi) - (chirp*1i-lambda0).*tfx0./tf0));
    
    for jcol = 1: N/2
        if abs(tf0(jcol)) > Threshold
            jcolhat = jcol - omega(jcol);
            
            if (jcolhat < Hidx + 1) && (jcolhat >= Lidx) &&...
                    (abs(jcolhat-Lidx+1-freq_idx(tidx)) < fband) %&& (abs(lambda(jcol)+ceil(N/2)-chirp_idx(tidx)) < cband)
                x_recon(tidx) = x_recon(tidx) + tf0(jcol);
            end
            
            
        end
    end
    
    
    
    
end

end
