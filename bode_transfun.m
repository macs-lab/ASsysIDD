function [magnitude,phase,f] = bode_transfun(tf,freq,SW_PLOT)
% function [magnitude,phase] = bode_transfun(tf,freq)
% computes the frequency response of a transfunction at indicated
% frequencies.
% input:    tf: transfer function
%           freq: frequency vector
% output:   magnitude
%           phase
if ~exist('freq')
    [compMag,compPhase,freq] = bode(tf);
else
    [compMag,compPhase] = bode(tf,freq*2*pi);
end
if ~exist('SW_PLOT','var')
    SW_PLOT = 1;
else
    SW_PLOT = 0;
end
f = freq;
magnitude = compMag(1,1,:);
magnitude = magnitude(:);
phase     = compPhase(1,1,:);
phase     = phase(:);
phase     = mod(phase+180,360)-180;
if SW_PLOT
    figure;
    subplot(211);
    plot(freq,20*log10(abs(magnitude)));
    grid on;
    ylabel('Gain (dB)');
    %      xlabel('Frequency (Hz)');
    %      axis('off')
    subplot(212);
    plot(freq,phase);
    grid on;
    ylabel('Phase (degree)');
    xlabel('Frequency (Hz)');
    
    % figure;
    % subplot(211);
    % semilogx(freq,20*log10(abs(magnitude)));
    % grid on;
    % ylabel('Gain (dB)');
    % %      xlabel('Frequency (Hz)');
    % %      axis('off')
    % subplot(212);
    % semilogx(freq,phase);
    % grid on;
    % ylabel('Phase (degree)');
    % xlabel('Frequency (Hz)');
end