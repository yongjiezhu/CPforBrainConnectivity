
%%%%  In order to run this code, please install the tensor toolboxes. It can be downloaded from
%%%%  https://www.tensortoolbox.org/
%%%%  Referenes for the toolboxes:
%%%%  Brett W. Bader, Tamara G. Kolda and others. MATLAB Tensor Toolbox,
%%%%  Version 3.1. Available online at https://www.tensortoolbox.org,
%%%%  2019

%%%%  This demo code shows how to use cp decompostion on time-frequency connnectivity data,
%%%%  and reproduce simulation in our TNSRE 2019 paper.
%%%%  References:
%%%%  Zhu, Y., Liu, J., Mathiak, K., Ristaniemi, T., & Cong, F. (2019). Deriving electrophysiological brain network connectivity via tensor component analysis during freely listening to music. IEEE Transactions on Neural Systems and Rehabilitation Engineering, 28(2), 409-418.
%%%%  Zhu, Y., Liu, J., Ye, C., Mathiak, K., Astikainen, P., Ristaniemi, T., & Cong, F. (2020). Discovering dynamic task-modulated functional networks with specific spectral modes using MEG. NeuroImage, 116924.
%%
addpath(genpath(fullfile(pwd,'tensor_toolbox')))
load simulatedData.mat
%%
R=3; % components
P = cp_als(tensor(double(Sim)),R); % cp decomposition using als algorithm
con=P.u{1}; % connectivity factor with dimension of paris of 68 nodes
spec=P.u{2}; % spectral factor
temporal=P.u{3}; % temporal factors
%% visualization
nNodes =68;
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
        %%% temporal factor
        code = (kcomp-1)*subplotNum+1;
        subplot(3,subplotNum,code)
        plot(tIndex,temporal(:,count),'linewidth',2)
        if kcomp==1
            title(['Temporal of feature #  ',int2str(count)],'fontsize',14)
            xlabel('Time/ms','fontsize',14)
            ylabel('Magnitude','fontsize',14)
        else
            title(['# ',int2str(count)],'fontsize',14)
        end
        %%% spectral factor
        code = (kcomp-1)*subplotNum+2;
        subplot(3,subplotNum,code)
        plot(fIndex,spec(:,count),'linewidth',2)
        grid on
        xlim([min(fIndex) max(fIndex)])
        if kcomp==1
            title(['Spectrum of feature #  ',int2str(count)],'fontsize',14)
            xlabel('Frequency/Hz','fontsize',14)
            ylabel('Magnitude','fontsize',14)
        else
            title(['# ',int2str(count)],'fontsize',14)
        end
        
        %%% Connectivity factor
        code = (kcomp-1)*subplotNum+3;
        subplot(3,subplotNum,code)
        AdjaMat=vect2conn(con(:,count),nNodes); % reshape connectivity vector to adjancent matrix
        imagesc(AdjaMat)
        if kcomp==1
            title(['Connectivity of feature # ',int2str(count)],'fontsize',13)
        else
            title(['# ',int2str(count)],'fontsize',13)
        end
    end
end
%%

