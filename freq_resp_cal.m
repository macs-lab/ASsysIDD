% This function computes the magnitude and phase of the frequency response between u (input)
% and y (output) using the spectral analysis.

function varargout=freq_resp_cal(y,u,Fs,Nfft)
% old output format:
%   function [mag,Freq]=freq_resp_cal(y,u,Fs)
% 2012-07-25: 
%   Xu Chen added the phase output

FRF = [];

%Fs = 800;
Ts = 1/Fs;
fmin_visu = 1;
fmax_visu = round(Fs/2);


if nargin < 4
Nfft = 2048/2;
else
end
window = hanning(Nfft);
noverlap = 0.75*Nfft;

%[Txy,Freq] = tfe(u,y,Nfft,Fs,window,noverlap) ;
[Txy,Freq] = tfestimate(u,y,window,noverlap,Nfft,Fs) ;



%
% figure;
% semilogx(Freq(ind),20*log10(abs(Txy(ind))));grid
% title('Magnitude of the frequency response of the secondary path obtained by spectral analysis')
% xlabel('Frequency (Hz)');
% ylabel('Magnitude (dB)');

% figure;
% plot(Freq(ind),180/pi*(angle(Txy(ind))));grid
% xlabel('Fréquency (Hz)');
% ylabel('Phase (°)');

mag=20*log10(abs(Txy));
pha=180/pi*(angle(Txy));
if nargout == 2
    varargout{1} = mag;
    varargout{2} = Freq;
elseif nargout == 3
    varargout{1} = mag;
    varargout{2} = Freq;
    varargout{3} = pha;
else   
    ind = find(Freq>fmin_visu & Freq<fmax_visu);
    figure;
    subplot(211)
    plot(Freq(ind),20*log10(abs(Txy(ind))));grid
    title('Frequency response of the secondary path obtained by spectral analysis')
    ylabel('Magnitude (dB)');
    subplot(212)
    plot(Freq(ind),pha(ind))
    grid
    xlabel('Frequency (Hz)');
    ylabel('Phase (deg)')
end