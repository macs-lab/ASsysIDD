% Readme:
%   required toolbox (should be added in your MATLAB path): 
%       landau's system id toolbox
%   This file conducts various recursive parameter adaptation algorithms to identify an active suspension plant
%   IMPORTANT: To adopt the file to your need, change the data set and orders of
%   tranfer functions, based on your specific physical hardware.
%
%   Make SURE that the mean values of u and y are zero first.
%
% Author: Xu Chen
load data_sec2;
% keyboard
if 0
    %%
    close all
end
u = u(1000:4048);
y = y(1000:4048);
TUNE_ORDER      = 0;
LOW_ORDER_ID    = 0;
if TUNE_ORDER
    if LOW_ORDER_ID
        START_ORDER = 16;
        END_ORDER   = 16;
    else
        START_ORDER = 10;
        END_ORDER   = 15;
    end
end
ALL_ID_METHOD   = 0;
TEST_ORDER      = 0;
figure,subplot(211),stairs(u),title('input')
subplot(212),stairs(y),title('output')
Ts = 1/800; Fs = 1/Ts;
[mag,freq] = freq_resp_cal(y,u,Fs);
if TEST_ORDER
    %% order test
    nmax = 26;
    % disp('======rls=====')
    % [V,S,VS] = estorderls(y,u,nmax)
    % figure,plot([1:nmax+1],V,[1:nmax+1],VS),legend('V','VS')
    disp('======iv======')
    [V,S,VS] = estorderiv(y,u,nmax)
    figure,plot([1:nmax+1],V,[1:nmax+1],VS),legend('V','VS')
else
    n = 22;
end
%%
if ALL_ID_METHOD
    %% rls
    [B_rls,A_rls] = rls(y,u,n,n,0)
    [wlossf,ulossf,wrni,urni,wyhat,uyhat]=olvalid(B_rls(1:end),A_rls,1,y,u);
    
    %% foloe
    % [B_foloe,A_foloe] = foloe(y,u,2,3,0,[1 -1 0.5])
    % [wlossf,ulossf,wrni,urni,wyhat,uyhat]=olvalid(B_foloe(1:end),A_foloe,1,y,u);
    %% oloe
    [B_oloe,A_oloe] = oloe(y,u,n,n,0)
    [wlossf,ulossf,wrni,urni,wyhat,uyhat]=olvalid(B_oloe(1:end),A_oloe,1,y,u);
end
%% afoloe
[B_afoloe,A_afoloe] = afoloe(y,u,n,n,0)
[wlossf,ulossf,wrni,urni,wyhat,uyhat]=olvalid(B_afoloe(1:end),A_afoloe,1,y,u);

tf_sec2_afoloe = tf(B_afoloe,A_afoloe,Ts);
[mag_afoloe,ph_afoloe,freq_afoloe] = bode_transfun(tf_sec2_afoloe,freq);
figure,pzplot(tf_sec2_afoloe)
%% xoloe
% n = 16
% finally chose 16th order regardless of the fact that sysID failed to pass
% the whitenness test. The fail of this test comes from the difference
% between the low freq systems. sec2 and sec3 are slightly different in
% this region
[B_xoloe,A_xoloe,C_xoloe] = xoloe(y,u,n,n,n,0,2000,1,1)
% [B_xoloe,A_xoloe,C_xoloe] = xoloe(y,u,n,n,n,0,1000,0.99,1)
% function [B,A,C]=xoloe(y,u,na,nb,nc,d,Fin,lam1,lam0)
[wlossf,ulossf,wrni,urni,wyhat,uyhat]=olvalid(B_xoloe(1:end),A_xoloe,C_xoloe,y,u);

tf_sec2 = tf(B_xoloe,A_xoloe,Ts);
[mag_xoloe,ph_xoloe,freq_xoloe] = bode_transfun(tf_sec2,freq);
figure,pzplot(tf_sec2)
B = B_xoloe;
A = A_xoloe;
% freq resp of the inverse plant
mag_invP = 1./mag_xoloe;
ph_invP = -ph_xoloe;
freq_invP = freq_xoloe;
if 0
    %%
    save data_sec2_blackbox_tf B A Ts
    save data_sec2_blackbox_inv_freq_resp  mag_invP ph_invP freq_invP
end
%% xoloe
if TUNE_ORDER
    for n = START_ORDER:END_ORDER
        % n = 10
        % finally chose 16th order regardless of the fact that sysID failed to pass
        % the whitenness test. The fail of this test comes from the difference
        % between the low freq systems. sec2 and sec3 are slightly different in
        % this region
        [B_xoloe_lowORD,A_xoloe_lowORD,C_xoloe_lowORD] = xoloe(y,u,n,n,n,0,2000,1,1)
        % [B_xoloe,A_xoloe,C_xoloe] = xoloe(y,u,n,n,n,0,1000,0.99,1)
        % function [B,A,C]=xoloe(y,u,na,nb,nc,d,Fin,lam1,lam0)
        [wlossf,ulossf,wrni,urni,wyhat,uyhat]=olvalid(B_xoloe_lowORD,A_xoloe_lowORD,C_xoloe_lowORD,y,u);
        
        tf_sec2_lowORD = tf(B_xoloe_lowORD,A_xoloe_lowORD,Ts);
        [mag_xoloe_lowORD,ph_xoloe_lowORD,freq_xoloe_lowORD] = bode_transfun(tf_sec2_lowORD,freq);
        figure,pzplot(tf_sec2_lowORD)
        title('low order approx')
        B = B_xoloe;
        A = A_xoloe;
        if 0
            %%
              save data_sec2_blackbox_tf_lowOrder B A Ts
        end
    end
end

%%
figure,hold on, grid on, zoom on
plot(freq,mag,...
    freq,20*log10(abs(mag_afoloe)),'g:.',...
    freq,20*log10(abs(mag_xoloe)),'r--');
title('Frequency responses of the measured and the identified systems')
legend('measured system','identified system: afoloe','identified system: xoloe')
%%
if LOW_ORDER_ID
    figure,hold on, grid on, zoom on
    plot(freq,mag,...
        freq,20*log10(abs(mag_xoloe)),'r--',...
        freq,20*log10(abs(mag_xoloe_lowORD)),'k:','Linewidth',2);
    title('Frequency responses of the measured and the identified systems')
    legend('measured system','identified system','system low order approximation')
else
    figure,hold on, grid on, zoom on
    plot(freq,mag,...
        freq,20*log10(abs(mag_xoloe)),'r--');
    title('Frequency responses of the measured and the identified systems')
    legend('measured system','identified system')
end