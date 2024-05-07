% Readme:
%   required toolbox (should be added in your MATLAB path): 
%       subspace_sysID toolbox shared in Google Drive
%   This file conducts subspace system identification algorithms to identify an active suspension plant
%   Read the files to understand the procedures. You will be promoted to
%   input the order during the execusion of the file.
%
%   Make SURE that the mean values of u and y are zero first.
%
% Author: Xu Chen
% 2012-07-25
% 
%
clear all
close all
load data_sec2
usave = u;
ysave = y;
% short summary 2012-07-25:
% 1, subspace sys ID seems to give more accurate results with the same
% system order. See the results on Landau's benchmark problem 2012
% 
% 2, if there are prior fixed parts in the transfer function
% e.g., double integrator or double differentiator type of systems
% it is suggested to try first compensate their effects through pre
% filtering. The reason is that if we rely the full order information all
% on the sys ID process, unstable zeros (especially near z=1 and
% z=-1) will mostly be identified to be outside of the unit circle. This
% causes problems for transfer-function inverse algorithms.
SW_SAVE = input('save?: yes-[1] no-[0 or Press Enter]: ');
if isempty(SW_SAVE)||SW_SAVE == 0
    SW_SAVE = 0;
else 
    SW_SAVE = 1;
end
N_TRUNK = 1500;
Ts = 1/800;
Fs = 1/Ts;
[mag,freq,pha] = freq_resp_cal(y,u,Fs,2000);
% [mag,freq,pha] = freq_resp_cal(y,u,Fs);
SW_INVERSE_ID = 0;
if SW_INVERSE_ID
    u = ysave(N_TRUNK:end);
    y = usave(N_TRUNK:end);
else
    SW_FIX_ZERO1 = 1; % prefilter the data to remove the fixed zero z=1
    if SW_FIX_ZERO1
        u = filter([1 -1],[1 0],u);
        u = filter([1 -1],[1 0],u);
    end
    u = u(N_TRUNK:end);
    y = y(N_TRUNK:end);
end

figure,subplot(211),stairs(u),title('input')
subplot(212),stairs(y),title('output')

%%
echo on;
max_ORDER = 40;
% max_ORDER = 26;
ORDER = 22;
min_ORDER = 16;
%   We will now identify this system from the data u,y
%   with the subspace identification algorithm: subid
%
%   The only extra information we need is the "number of block rows" i
%   in the block Hankel matrices.  This number is easily determined
%   as follows:

%   Say we don't know the order, but think it is maximally equal to 26.
%
%       max_order = 26;
%
%   As described in the help of subid we can determine "i" as follows:
%
%       i = 2*(max_order)/(number of outputs)
%
i = 2*(max_ORDER)/1;

%   Hit any key
% pause
clc
%
%   The subspace algorithms is now easily started.
%   When prompted for the system order you should enter 22.
%
[A,B,C,D] = subid(y,u,i);

%   Did you notice the order was very easy to determine
%   from the singular values?

%   Hit any key
% pause
clc

w = [0:0.005:0.5]*(2*pi); 		% Frequency vector
%   Just to make sure we identified the original system again,
%   we will compare the original and estimated transfer function.
%
if 0
    M1 = dbode(A,B,C,D,1,1,w);
    %     M2 = dbode(A,B,C,D,1,2,w);
    figure(1)
    hold off;
    %     subplot;clg;
    %     subplot(221);plot(w/(2*pi),[m1(:,1),M1(:,1)]);title('Input 1 -> Output 1');
    plot(w/2/pi,M1(:,1))
    %     subplot(222);plot(w/(2*pi),[m2(:,1),M2(:,1)]);title('Input 2 -> Output 1');
    %     subplot(223);plot(w/(2*pi),[m1(:,2),M1(:,2)]);title('Input 1 -> Output 2');
    %     subplot(224);plot(w/(2*pi),[m2(:,2),M2(:,2)]);title('Input 2 -> Output 2');
    
    %   As you can see, the original and identified system are
    %   exactly the same.
    %
    %   Hit any key
    pause
    clc
end

%%
PsubID = ss(A,B,C,D,Ts);
if SW_FIX_ZERO1
    PsubID = PsubID*tf([1 -2 1],[1 0 0],Ts);
end
ssModel_sec = PsubID;
if SW_SAVE
    save landau_ssModel_sec ssModel_sec
end
[mag_sub,ph_sub,freq_sub] = ...
    bode_transfun(PsubID,freq,0);

figure, pzplot(PsubID)
figure,
subplot(211)
hold on, grid on, zoom on
plot(freq,mag,...
    freq,20*log10(abs(mag_sub)),'r--');
title('Frequency responses of the measured and the identified systems')
subplot(212)
plot(freq,pha,...
    freq,ph_sub,'r--');
legend('measured system','identified system')
%%
%%
%   The function "simul" allows you to check the size of the simulation
%   error.  this is a measure for the difference between the original
%   and the simulated output:
%
[ys,ers] = simul(y,u,A,B,C,D);
%
%   ers contains the error per output in percentage:
ers
%   While ys contains the simulated output:
figure
plot([y(100:400,1),ys(100:400,1)])
title('Real (yellow) and simulated (purple) output')

%   They coincide well.
%%
if SW_INVERSE_ID
    return
end
%%
Pl = importdata('landau_model_sec.mat');
% Pl = importdata('new_model_sec.mat');
Plandau = tf(Pl.B,Pl.A,Pl.Ts,'variable','z^-1');
[mag_giv,ph_giv,freq_giv] = ...
    bode_transfun(Plandau,freq,0);
figure,
subplot(211)
hold on, grid on, zoom on
plot(freq,mag,...
    freq,20*log10(abs(mag_sub)),'r--',...
    freq,20*log10(abs(mag_giv)),'g:');
ylim([-70,30])
title('Frequency responses of the measured and the identified systems')
legend('measured system','subspace identified system',...
    'recursive identified system')
subplot(212)
plot(freq,pha,...
    freq,ph_sub,'r--',...
    freq,(ph_giv),'g:');
%%
return
%%
%   Hit any key
pause
clc
%
%   The subspace identification function subid also allows you to
%   identify the noise system.  Note that the fourth argument in the
%   following subid call is the system order.  When this is given, the
%   singular values are not plotted and you are not prompted for the order.

[A,B,C,D,K,R] = subid(y,u,i,ORDER);

%   Hit any key
pause
clc
%
%   you can now compute the "one step ahead" prediction using predic:

[yp,erp] = predic(y,u,A,B,C,D,K);

%   Compare the prediction error and simulation error:
[ers;erp]

%   The prediction errors are significantly smaller.

%   Hit any key
pause
clc

%   In many practical examples, the gap in the singular value plot
%   is not as clear as in this example.  The order decision then becomes
%   less trivial.  There is however an nice feature of subid which allows
%   for fast computation of systems with different orders.
%   This can be done through the extra variable AUX which appears
%   as an input as well as an output argument of subid.
%   The last parameter (1) indicates that the algorithm should run silently.
%
[A,B,C,D,K,R,AUX] = subid(y,u,i,2,[],[],1);
era = [];
for n = min_ORDER:max_ORDER
    [A,B,C,D,K,R] = subid(y,u,i,n,AUX,[],1);
    [ys,ers] = simul(y,u,A,B,C,D);
    era(n-min_ORDER+1,:) = ers;
end

%   Hit any key
pause
clc
%
%   We have now determined the simulation errors for all systems
%   from order 1 through 6.
%   Plotting these errors often gives a clearer indication of the order:
%
figure;
bar([min_ORDER:max_ORDER],era);
%     bar(era)
title('Simulation error');
xlabel('System order');

%   It now even becomes more clear that the system order is 4.
%
%   Hit any key
pause
clc
%
%   We did find this very useful, so we included the above code in
%   a function: allord.  The above result could have been obtained
%   with the one line code:
%
[ersa,erpa] = allord(y,u,i,[min_ORDER:max_ORDER],AUX);

%   Hit any key
pause
clc

%   A last feature we would like to illustrate is that subid also
%   can work with principal angles instead of singular values.
%   The order of the system is then equal to the number of principal
%   angles different from 90 degrees:

[A,B,C,D,K,R] = subid(y,u,i,[],AUX,'cva');

%   The order is still clearly equal to 4.
%
%   This concludes the startup demo.  You should now be ready
%   to use subid on your own data.
%
%   Note that you can also identify time series (no input) with subid.
%   See the help of subid and sto_demo for more explanation.

echo off