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
