%% LOADING DATA 
f_samp= 2; %Hz
h = cell(1,5); % Contains Healthy HRVs
p = cell(1,5); % Contains IUGS Patological HRVs
for i=1:5
    load(['Healthy_' num2str(i)]);
    data = 1./(data./60)*1000; % Conversion from bpm to milliseconds
    data(data==Inf) = NaN;     % Setting "red values" to NaN
    [~,TF] = rmoutliers(data); % Removing outliers by setting them to NaN
    data(TF) = NaN;            %                    "
    h{i} = data;
    load(['IUGR_' num2str(i)]);
    data = 1./(data./60)*1000;
    data(data==Inf) = NaN;
    [~,TF] = rmoutliers(data);
    data(TF) = NaN;
    p{i} = data;
    clear data
end

%% PREPROCESSING
h_flat = cell(1,5);
p_flat = cell(1,5);

movmeansamps = 50;
% Baseline removal
for i=1:5
   h_flat{i} = h{i} - mean(h{i},'omitnan');
   p_flat{i} = p{i} - mean(p{i},'omitnan');
end

h_clean = cell(1,5);
p_clean = cell(1,5);
nwindows = 50;      % Number of windows chosen for data interpolation
nan_threshold = 10; % Number of consecutive NaNs for the window to be removed

% Signal truncation (up to a number of samples multiple of 'nwindows')
for i=1:5
    h_flat{i} = h_flat{i}(1:floor(length(h_flat{i})/nwindows)*nwindows);
    p_flat{i} = p_flat{i}(1:floor(length(p_flat{i})/nwindows)*nwindows);
end

% Signal cleaning from NaN values and interpolation
for i=1:5
   h_clean{i}= clean_nans(h_flat{i}, nwindows, nan_threshold); % See function
   p_clean{i}= clean_nans(p_flat{i}, nwindows, nan_threshold); % See function
   h_clean{i} = fillmissing(h_clean{i},'linear'); %Interpolation of NaN values
   p_clean{i} = fillmissing(p_clean{i},'linear'); %            "
end


%% SPECTRAL ANALYSIS
nfft = 1024; % Number of samples used for FFT
h_Bartlett = cell(1,5);
p_Bartlett = cell(1,5);
h_Welch = cell(1,5);
p_Welch = cell(1,5);
h_YW = cell(1,5);
p_YW = cell(1,5);
for i=1:5
    % Bartlett's PSD using implicit rectangular window of double length
    % than the one chosen for the preprocessing
    winwidth = length(h_clean{i})/nwindows*2;
    [h_Bartlett{i}, fspace] = pwelch(h_clean{i},rectwin(winwidth),0,nfft,f_samp);
    [p_Bartlett{i}, ~] = pwelch(p_clean{i},rectwin(winwidth),0,nfft,f_samp);
    
    overlap = 0.5; % Percentage of window overlap
    % Welch PSD using explicit hamming windows and 50% overlap
    [h_Welch{i}, ~] = pwelch(h_clean{i},hamming(winwidth),winwidth*overlap,nfft,f_samp);
    [p_Welch{i}, ~] = pwelch(p_clean{i},hamming(winwidth),winwidth*overlap,nfft,f_samp);
    
    % Yule-Walker Parametric PSD with order
    order = 12  ;
    [h_YW{i}, ~] = pyulear(h_clean{i},order,nfft,f_samp);
    [p_YW{i}, ~] = pyulear(p_clean{i},order,nfft,f_samp);
end

%% SPECTRAL DECOMPOSITION
%Definition of frequency intervals
VLF = fspace<0.03;
LF = fspace>=0.03 & fspace<0.15;
MF = fspace>=0.15 & fspace<0.5;
HF = fspace>=0.5;
hpow = zeros(5,4);
ppow = zeros(5,4);

% Spectral decomposition of every signal with Trapezolidal Mehtod
for i=1:5
    hpow(i,1) = trapz(h_Welch{i}(VLF));
    hpow(i,2) = trapz(h_Welch{i}(LF));
    hpow(i,3) = trapz(h_Welch{i}(MF));
    hpow(i,4) = trapz(h_Welch{i}(HF));
    ppow(i,1) = trapz(p_Welch{i}(VLF));
    ppow(i,2) = trapz(p_Welch{i}(LF));
    ppow(i,3) = trapz(p_Welch{i}(MF));
    ppow(i,4) = trapz(p_Welch{i}(HF));
end

%% PARAMETERS CALCULATIONS
%Calculating the ratio = LF/(MF+HF)
h_ratio = hpow(:,2)./(hpow(:,3)+hpow(:,4));
p_ratio = ppow(:,2)./(ppow(:,3)+ppow(:,4));
h_apen = zeros(5,1);
p_apen = zeros(5,1);

% APEN takes a lot of time. If available, load the ones altready calculated
load('entropies.mat');

% for i=1:5
%     % Calculating ApEn
%     m = 3; % Standard value across literature
%     r_h = 0.15*std(h_clean{i}); % r = 0.15*std(signal)
%     r_p = 0.15*std(p_clean{i});
%     h_apen(i) = ApEn(h_clean{i},m,r_h); % See function
%     p_apen(i) = ApEn(p_clean{i},m,r_p);
% end

t = table([h_ratio; p_ratio], [h_apen; p_apen], 'VariableNames', {'LFHF_ratio','ApEn'});

% All results are saved in a struct
results.LFHF_mean_healthy = mean(t.LFHF_ratio(1:5));
results.LFHF_mean_patho = mean(t.LFHF_ratio(6:end));
%% CSA CALCULATION
    hpsdwind = cell(1,5);
    ppsdwind = cell(1,5);
    for i=1:5
        w = 1;
        for j=1:(length(h_clean{i})/nwindows):length(h_clean{i}) % Scans each windows
            K = j:j+length(h_clean{i})/nwindows-1;
            %PSD of each window
            hpsdwind{i}(w,:) = fft(h_clean{i}(K),nfft).*conj(fft(h_clean{i}(K),nfft))/length(K)';
            w=w+1;
        end
    end
    for i=1:5
        w = 1;
        for j=1:(length(p_clean{i})/nwindows):length(p_clean{i}) % Scans each windows
            K = j:j+length(p_clean{i})/nwindows-1;
            ppsdwind{i}(w,:) = fft(p_clean{i}(K),nfft).*conj(fft(p_clean{i}(K),nfft))/length(K)';
            w=w+1;
        end
    end

%% ONE-FIGURE PLOT
close all
for i=1:5
    t = (1:length(h{i}))/2;
    figure('Name',['Healthy ' num2str(i)],'units','normalized',...
        'outerposition',[0.05 0.05 0.9 0.9])
    subplot('position',[0.0553,0.7413,0.386,0.1836]); % Raw signal
    plot(t,h{i}); title(['Healthy ' num2str(i) ' - RAW Signal']);
    xlabel('Time [sec]');
    ylabel('bpm [adim.]'); 
%     set(gca,'Position',[0.0553,0.7413,0.386,0.1836]);
    t = (1:length(h_clean{i}))/2;
    subplot('position',[0.0553,0.450059453032105,0.386,0.1836]); % Cleaned signal
    plot(h_clean{i}); title(['Healthy ' num2str(i) ' - CLEAN Signal']);
    xlabel('Time [sec]');
%     set(gca,'Position',[0.065392731535756,0.450059453032105,0.386869871043376,0.183620689655173]);
    
    subplot('position',[0.0553,0.11,0.386,0.2544]); % PSD - Power Spectral Density
    plot(fspace, h_Welch{i}); 
    hold on; plot(fspace, h_Bartlett{i}); plot(fspace, h_YW{i}); 
    title(['Healthy ' num2str(i) ' - Power Spectral Density PSD']);
    legend('Welch','Bartlett','Yule-Walker');
    xlabel('Frequency [Hz]');
    ylabel('Magnitude [Absolute Units]');
%     set(gca,'Position',[0.055392731535756,0.11,0.387456037514654,0.254447086801427]);
    
    subplot('position',[0.501778022698477,0.11,0.439076577791433,0.81143318094623]); % CSA - Compressed Spectral Analysis
%     set(gca,'Position',[0.501778022698477,0.11,0.439076577791433,0.81143318094623]);
    [F,N] = meshgrid(1:size(hpsdwind{i},2),size(hpsdwind{i},1):-1:1);
    F = F/nfft*f_samp;
    mesh(N,F,hpsdwind{i});
    ylim([0 F(end/2)]);
    ylabel('Frequency [Hz]');
    xlabel('Window [k]');
    xticklabels(num2cell(string(flip(xticks))))
    zlabel('PSD Magnitude [Absolute Units]');
    title(['Healthy ' num2str(i) ' - Compressed Spectral Analysis CSA']);
    xlim([0 nwindows]);
    ca = gca;
    ca.Children.MeshStyle = 'row';
    ca.Children.EdgeColor = 'interp';
    ca.Children.FaceColor = 'interp';
    ca.Children.FaceAlpha = 0.2;
    ca.Children.LineWidth = 1;
    view([119 52]);
    colormap jet
    
%     subplot('position',[0.781260504201681,0.3,0.123739495798319,0.815]); % Piechart
%     pie(hpow(i,:)); title('Spectral Decomposition');
end

for i=1:5
    t = (1:length(p{i}))/2;
    
    figure('Name',['IUGR ' num2str(i)],'units','normalized',...
        'outerposition',[0.05 0.05 0.9 0.9])
    subplot('position',[0.0553,0.7413,0.386,0.1836]); % Raw signal
    plot(t,p{i}); title(['IUGR ' num2str(i) ' - RAW Signal']);
        xlabel('Time [sec]');
    ylabel('bpm [adim.]'); 
%     set(gca,'Position',[0.0553,0.7413,0.386,0.1836]);
    t = (1:length(p_clean{i}))/2;
    subplot('position',[0.0553,0.450059453032105,0.386,0.1836]); % Cleaned signal
    plot(t, p_clean{i}); title(['IUGR ' num2str(i) ' - CLEAN Signal']);
        xlabel('Time [sec]');
%     set(gca,'Position',[0.065392731535756,0.450059453032105,0.386869871043376,0.183620689655173]);
    
    subplot('position',[0.0553,0.11,0.386,0.2544]); % PSD - Power Spectral Density
    plot(fspace, p_Welch{i}); 
    hold on; plot(fspace,p_Bartlett{i}); plot(fspace, p_YW{i}); 
    title(['IUGR ' num2str(i) ' - Power Spectral Density PSD']);
    legend('Welch','Bartlett','Yule-Walker');
%     set(gca,'Position',[0.055392731535756,0.11,0.387456037514654,0.254447086801427]);
    
    subplot('position',[0.501778022698477,0.11,0.439076577791433,0.81143318094623]); % CSA - Compressed Spectral Analysis
%     set(gca,'Position',[0.501778022698477,0.11,0.439076577791433,0.81143318094623]);
    [F,N] = meshgrid(1:size(ppsdwind{i},2),size(ppsdwind{i},1):-1:1);
    F = F/nfft*f_samp;
    mesh(N,F,ppsdwind{i});
    ylim([0 F(end/2)]);
    ylabel('Frequency [Hz]');
    xlabel('Window [k]');
    xticklabels(num2cell(string(flip(xticks))))
    zlabel('PSD Magnitude [Absolute Units]');
    title(['IUGR ' num2str(i) ' - Compressed Spectral Analysis CSA']);
    xlim([0 nwindows]);
    ca = gca;
    ca.Children.MeshStyle = 'row';
    ca.Children.EdgeColor = 'interp';
    ca.Children.FaceColor = 'interp';
    ca.Children.FaceAlpha = 0.2;
    ca.Children.LineWidth = 1;
    view([119 52]);
    colormap jet
    
%      subplot('position',[0.781260504201681,0.3,0.123739495798319,0.815]); % Piechart
%      pie(ppow(i,:)); title('Spectral Decomposition');
end

%% PIE CHARTS
%Healthy
hhfigure
for i=1:5
    subplot(2,5,i); pie(hpow(i,:)); title(['Healthy ' num2str(i)]);
    legend
end
for i=1:5
    subplot(2,5,i+5); pie(ppow(i,:)); title(['Patho ' num2str(i)]);
end

%% SIGNAL PLOTTING
% Healthy
for i=1:5
    figure
    subplot(3,1,1);
    plot(h{i});
    xticklabels(xticks/2);
    xlabel('Time [sec]');
    title(['GIVEN signal - Healthy ' num2str(i)]);
    grid on
    subplot(3,1,2);
    plot(h_clean{i});
    xticklabels(xticks/2);
    xlabel('Time [sec]');
    title(['CLEANED signal - Healthy ' num2str(i)]);
    grid on
    subplot(3,1,3);
    plot(fspace,h_Bartlett{i});
    hold on
    plot(fspace,h_Welch{i});
    plot(fspace,h_YW{i});
    xlabel('Frequency [Hz]');
    title(['PSD - Healthy ' num2str(i)]);
    grid on
    legend('Bartlett',['Welch [Hamming ' num2str(overlap*100) '% overlap'],['Yule-Walker [Ord. ' num2str(order) ']']);
end
%Pathological
for i=1:5
    figure
    subplot(3,1,1);
    plot(p{i});
    xticklabels(xticks/2);
    xlabel('Time [sec]');
    title(['GIVEN signal - Pathological ' num2str(i)]);
    grid on
    subplot(3,1,2);
    plot(p_clean{i});
    xticklabels(xticks/2);
    xlabel('Time [sec]');
    title(['CLEANED signal - Pathological ' num2str(i)]);
    grid on
    subplot(3,1,3);
end
%% CSA PLOT
hpsdwind = cell(1,5);
ppsdwind = cell(1,5);
for i=1:5
    w = 1;
    for j=1:(length(h_clean{i})/nwindows):length(h_clean{i}) % Scans each windows
        K = j:j+length(h_clean{i})/nwindows-1;
        %PSD of each window
        hpsdwind{i}(w,:) = fft(h_clean{i}(K),nfft).*conj(fft(h_clean{i}(K),nfft))/length(K)';
        w=w+1;
    end
 %Plotting the CSA 
    figure
    [F,N] = meshgrid(1:size(hpsdwind{i},2),1:size(hpsdwind{i},1));
    F = F/nfft*f_samp;
    mesh(N,F,hpsdwind{i});
    ylim([0 F(end/2)]);
    ylabel('Frequency [Hz]');
    xlabel('Window [k]');
    zlabel('PSD');
    title(['Healthy ' num2str(i) ' CSA']);
    xlim([0 nwindows]);
    ca = gca;
    ca.Children.MeshStyle = 'row';
    ca.Children.EdgeColor = 'interp';
    ca.Children.FaceColor = 'interp';
    ca.Children.FaceAlpha = 0.2;
    ca.Children.LineWidth = 1;
    view([119 52]);
    colormap jet
    
    w = 1;
    for j=1:(length(p_clean{i})/nwindows):length(p_clean{i}) % Scans each windows
        K = j:j+length(p_clean{i})/nwindows-1;
        ppsdwind{i}(w,:) = fft(p_clean{i}(K),nfft).*conj(fft(p_clean{i}(K),nfft))/length(K)';
        w=w+1;
    end
    figure
    [F,N] = meshgrid(1:size(ppsdwind{i},2),1:size(ppsdwind{i},1));
    F = F/nfft*f_samp;
    mesh(N,F,ppsdwind{i});
    ylim([0 F(end/2)]);
    ylabel('Frequency [Hz]');
    xlabel('Window [k]');
    zlabel('PSD');
    title(['Pathological ' num2str(i) ' CSA']);
    xlim([0 nwindows]);
    ca = gca;
    ca.Children.MeshStyle = 'row';
    ca.Children.EdgeColor = 'interp';
    ca.Children.FaceColor = 'interp';
    ca.Children.FaceAlpha = 0.2;
    ca.Children.LineWidth = 1;
    view([119 52]);
    colormap jet
end

%% HEALTHY/IUGR WELCH PSD COMPARISON
figure; hold on; grid on
plot(fspace,h_Welch{1},'Color','#0000ff')
plot(fspace,h_Welch{2},'Color','#729bf2')
plot(fspace,h_Welch{3},'Color','#72d5f2')
plot(fspace,h_Welch{4},'Color','#6980f0')
plot(fspace,h_Welch{5},'Color','#4fc3c0')
plot(fspace,p_Welch{1},'Color','#c3554f')
plot(fspace,p_Welch{2},'Color','#e55202')
plot(fspace,p_Welch{3},'Color','#ff0000')
plot(fspace,p_Welch{4},'Color','#c58202')
plot(fspace,p_Welch{5},'Color','#dd9772')
xlim([0 0.5]);
title('Welch PSD comparison for Healthy and IUGR');
legend('Healthy 1','Healthy 2','Healthy 3','Healthy 4','Healthy 5',...
    'IUGR 1','IUGR 2','IUGR 3','IUGR 4','IUGR 5');

%% ENTROPY COMPARISON
figure; bar([t.ApEn(1:5)' ; t.ApEn(6:10)']')
title('ApEn calculation');

%% FUNCTIONS IMPLEMENTED
function cleaned_sig= clean_nans(sig, nwindows, thr)
%CLEAN_NANS(sig,nwindows,thr) outputs the signal cleaned from intervals
%with an amount of NaNs specified by 'thr'. The interval depends on the
%choice of the number of windows 'nwindows'
    for j=1:(length(sig)/nwindows):length(sig) % Scans each windows
        K = j:j+length(sig)/nwindows-1;
        for k=j:j+length(sig)/nwindows-1 % Scans each element in the window
            e=0; % Number of consecutive NaN elements
            while isnan(sig(k+e))
                e = e+1;
                if (k+e) > length(sig) 
                    break; 
                end
            end
            if e >= thr
                sig(K) = 0;
                break;
            end
        end
    end
    cleaned_sig = sig;
end

function apen = ApEn(x,m,r)
%APEN(x,m,r) calculates the Approximate Entropy (Pincus, 1990) associated
%to the signal 'x' and evaluated on 'm' samples with threshold 'r'.
    N = length(x);
    p = zeros(N-m+1,m);
    for i=1:N-m+1
        p(i,:) = x(i:i+m-1);
    end
    c = zeros(1,N-m+1);
    for i=1:N-m+1
        for j=1:N-m+1
            c(i) = c(i)+all(abs(p(i,:)-p(j,:))<r);
        end
    end
    c = c/(N-m+1);
    apen = sum(log(c))/(N-m+1);
    m=m+1;
    p = zeros(N-m+1,m);
    for i=1:N-m+1
        p(i,:) = x(i:i+m-1);
    end
    c = zeros(1,N-m+1);
    for i=1:N-m+1
        for j=1:N-m+1
            c(i) = c(i)+all(abs(p(i,:)-p(j,:))<r);
        end
    end
    c = c/(N-m+1);
    apen = apen-sum(log(c))/(N-m+1);
end

