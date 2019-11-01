
% Author: Yongjie Zhu and Fengyu Cong
% Email: yongjie.zhu@foxmail.com
% Affiliation: Dalian University of Technology, China
%              University of Jyväskylä, Finland
% Date: October 27th, 2019.
%
% Require MATLAB Tensor Toolbox from
% http://www.sandia.gov/~tgkolda/TensorToolbox/
%
%% demo for simulated data
load simulatedData
rankR = 3; % the number of compnents

% constructing tensor
A{1} = Connection;  % Connectivity factors(Motor, Visual,pari-frontal)
A{2} = Spectral;    % Spectral factors (3,8,15 Hz)
A{3} = T;           % Temporal factors (sawtooth,square,sin)
Y = ktensor(A(:));  % a Kruskal tensor of rank-R
orgY = Y;

% Add Gaussian noise
SNR = 0; % Noise level in dB, inf for noiseless tensor
normY = norm(Y);
if ~isinf(SNR)
    Y = full(Y);
    sig2 = normY./10.^(SNR/20)/sqrt(prod(size(Y)));
    Y = Y + sig2 * randn(size(Y));
end
Y = tensor(Y); % For the real data, tensor Y can be obtained by
% computing the time-frequency phase-coupling (e.g.,wPLI,PLI)
% or envelope correlation based on sliding window
% between pairs brain regions.
%% run CP decomposition

fprintf(['Decomping' '\n'])
options = struct('tol',1e-8,'maxiters',500,'init','random',...
    'orthoforce',0,'fitmax',.99999,'verbose',1); %%%% set tol maxiters
R=3;
[P,Ah,fitarr] = ntf_fastHALS(Y,R,options); % can also be ALS()in Tensor Toolbox
con=P.u{1};
spec=P.u{2};
temporal=P.u{3};
%% some plots
subplotNum = 3;
figNum = ceil(R/3);
count = 0;
for is = 1:figNum
    figure
    set(gcf,'outerposition',get(0,'screensize'))
    for kcomp = 1:subplotNum
        count = count + 1;
        if count>R
            break;
        end
        %%% temporal component
        code = (kcomp-1)*subplotNum+1;
        subplot(3,subplotNum,code)
        plot(tIndex,temporal(:,count),'linewidth',2)
        if kcomp==1
            title(['Temporal profile #  ',int2str(count)],'fontsize',14)
            xlabel('Time/S','fontsize',14)
            ylabel('Magnitude','fontsize',14)
        else
            title(['# ',int2str(count)],'fontsize',14)
        end
        %%% spectral component
        code = (kcomp-1)*subplotNum+2;
        subplot(3,subplotNum,code)
        plot(fIndex,spec(:,count),'linewidth',2)
        grid on
        xlim([min(fIndex) max(fIndex)])
        if kcomp==1
            title(['Spectral profile #  ',int2str(count)],'fontsize',14)
            xlabel('Frequency/Hz','fontsize',14)
            ylabel('Magnitude','fontsize',14)
        else
            title(['# ',int2str(count)],'fontsize',14)
        end
        
        %%% spatial component
        code = (kcomp-1)*subplotNum+3;
        subplot(3,subplotNum,code)
        EdgeSize=vect2conn(con(:,count),68);
        imagesc(EdgeSize);
        if kcomp==1
            title(['Connectivity profile # ',int2str(count)],'fontsize',13)
            xlabel('Node#')
            ylabel('Node#')
        else
            title(['# ',int2str(count)],'fontsize',13)
        end
    end
end