%**************************************************************************
% Noise bypass and event detection from RAW data
% Author: Gabriel Moreira
% Reviewed: 14 June 2018
%**************************************************************************

%% Initialization
clc
clear all;

%% Noise Suppresion

% Loads image
load kernel.mat;
integration_time = 120000;
img_name = ['image_ms_' num2str(integration_time) '.jpg'];
frame_raw = imread(img_name);
jpgraw;
[height_raw, width_raw] = size(frame_raw);

%% Plots raw frame without any treatment

figure('pos',[10 10 900 600])
surf(frame_raw(1:200,1:200));
set(gca,'fontsize',20)
xlabel('x axis (pixel)', 'fontsize',20);
ylabel('y axis (pixel)', 'fontsize', 20);
zlabel('ADC count / Ionization charge (A.U.)', 'fontsize', 20);
view([45,45]);
title(['Raw subframe (integrated for ', num2str(integration_time), ...
    'ms)'], 'fontsize',20);
axis([0 200 0 200 0 1050]);

%% Applies wiener filter for noise reduction

adptv_noise_correction = wiener2(frame_raw,[10 10]);
figure('pos',[10 10 900 600])
surf(adptv_noise_correction(1:800,1:800));
view([45,45]);
title('Subframe after applying a Wiener filter', 'fontsize',20);
xlabel('x axis (pixel)', 'fontsize',20);
ylabel('y axis (pixel)', 'fontsize', 20);
axis([0 200 0 200 0 1050]);

%% Applies a filter whose kernel is an event

kernel = kernel / mean(mean(kernel));
adptv_noise_correction = imfilter(adptv_noise_correction, kernel);
figure('pos',[10 10 900 600]);
surf(adptv_noise_correction(1:200,1:200));
view([45,45]);
title('Subframe after convolution with a kernel', 'fontsize',20);
xlabel('x axis (pixel)', 'fontsize',20);
ylabel('y axis (pixel)', 'fontsize', 20);
axis([0 200 0 200 0 30000]);

%% Conditional Probability
% Probability of a pixel having a value higher than the threshold t,
% given that at least one adjacent pixel verifies this condition

v = adptv_noise_correction;
t = 80;
[ii, jj] = find(frame_raw >= t);
prob = length(ii)/(width_raw*height_raw)
cumsum = 0;

for p = 1:length(ii)
        if (ii(p) + 1 <= height_raw && jj(p) + 1 <= width_raw && ...
                ii(p) - 1 >= 1 && jj(p) - 1 >= 1)
            if (v(ii(p)+1, jj(p)) >= t || v(ii(p)-1, jj(p)) >= t || ...
                    v(ii(p), jj(p)+1) >= t || v(ii(p), jj(p)-1) >= t)
                cumsum = cumsum + 1;
            end
        end
end

condprob = cumsum/length(ii)

%% Computes new matrix with the standard deviations of the N*M submatrices

N = 8;
M = 8;
std_frame = SubDeviations(N, M, adptv_noise_correction);
std_dev_mean = mean(mean(std_frame));

% This is the frame that helps the identification of particles
id_frame = (std_frame / std_dev_mean);

%% Determines which submatrices have spots (could be events or not)

centers = [];
[h, w] = size(id_frame);
sorted_id = sort(id_frame(:), 'descend');
ratio_threshold = 0.4*mean(sorted_id(1:30));

for i = 1:h
    for j = 1:w
        if id_frame(i,j) > ratio_threshold
            centers = [centers; i*N j*M];
        end
    end
end

%% Retrieves a noise sample from the raw frame, knowing the spots positions

copy = frame_raw;
noise_sample = [];
k = 1;

for c = 1:length(centers)
    X = centers(c, 1);
    Y = centers(c, 2);
    copy(X:X+N, Y:Y+M) = 0;
end
for i = floor(1 + height_raw*rand(1,1000)) 
    for j = floor(1 + width_raw*rand(1,1000))
        if copy(i,j) ~= 0
            noise_sample(k) = copy(i, j);
            k = k + 1;
        end
    end
end

clearvars copy;
noise_mean = mean(noise_sample);
noise_stdvar = sqrt(var(noise_sample));
correction_frame = normrnd(-noise_mean, noise_stdvar, size(frame_raw));
frame_raw = frame_raw + correction_frame;

%% Calculates correlation between adjacent cells 

cumsum = 0;
for i = 2:801
     for j = 2:801
        cumsum = (1/800^2)*(frame_raw(i,j) - noise_mean)*(frame_raw(i-1,j)...
            -  noise_mean);
     end
end
corr = cumsum/(noise_stdvar^2);

%% Determines if the spot is an event. If so, it stores them in "events"

n_spots = length(centers);
event_size = 7; %-- event's maximum size
n_events = 0;
maximum = zeros(1, n_spots);
max_list = [];
events = zeros(event_size, event_size);
raw_x = []; %-- array of positions in the raw frame
raw_y = []; %-- array of positions in the raw frame

% Generate a copy of frame_raw with borders of thickness N and M
A = frame_raw;
A = [zeros(height_raw, N) A zeros(height_raw, N)];
[z, d] = size(A);
A = [zeros(M, d); A; zeros(M, d)];
ADU_threshold = 300; %-- adjust accordingly

for s = 1:n_spots
    spot = [];
    i = centers(s, 1);
    j = centers(s, 2);
    ilower = i;
    iupper = i+N-1;
    jlower = j;
    jupper = j+M-1;
    if iupper > height_raw
        iupper = height_raw;
    end
    if jupper > width_raw
        jupper = width_raw;
    end
    spot = frame_raw(ilower:iupper,  jlower:jupper);
    maximum(s) = max(max(spot));  
    [X,Y] = find(spot==maximum(s));
    
    % redefine the spot in the frame_raw coordinates
    raw_x(s) = X + ilower;
    raw_y(s) = Y + jlower;
    
    % spot centered in the pixel with max reading
    spot = A((raw_x(s)-floor(event_size/2)+M-1):(raw_x(s)+ ...
        floor(event_size/2)+M-1),(N+raw_y(s)-floor(event_size/2)-1): ...
        (N+raw_y(s)+floor(event_size/2))-1);
    
    % Get the neighbouring pixels and sort them
    surrounding = spot(floor(event_size/2):floor(event_size/2)+2, ...
        floor(event_size/2):floor(event_size/2)+2);
    sort_surrounding = sort(surrounding(:), 'descend');
    
    % Check if spot is an event:
    % Must have max higher than threshold
    % Must have at least a neighbouring pixel higher than threshold
    % Must be at least 3px away from nearest maximum
    if maximum(s) > ADU_threshold && sort_surrounding(2) > ADU_threshold
        
        for w = 1:length(n_events)
            if norm([raw_x(s), raw_x(s)] - [raw_x(w), raw_x(w)]) < 3 ...
                    && (i ~= w)
                break;
            end
            events(:,:,n_events + 1) = spot;
            n_events = n_events + 1;
            max_list(n_events + 1) = maximum(s);
        end
        
    end
end

%% Energy for each spot and average electron count

CVF = 0.2351; %-- check if true;
electrons = zeros(1, n_events);
electron_hole_si = 3.6; %-- eV
for c = 1:n_events
    ADU_sum = 0;
    for i = 1:event_size 
        for j = 1:event_size
            if events(i,j,c) > 0
                ADU_sum = ADU_sum + uint64(events(i,j,c));
            end
        end
    end
    electrons(c) = ADU_sum/CVF; %-- eV
end
avg_electrons = mean(electrons);
total_electrons = sum(electrons);
disp(['>> AVERAGE N. ELECTRONS PER EVENT: ' ...
    num2str(avg_electrons) ' eV']);

%% DATA

data.events = events;
data.n_events = n_events;
data.max = max_list;
data.avg_count = avg_electrons;
data.mean_noise = noise_mean;

%% Plot event's surfaces and saves the figures as PNG

fig = figure('units','normalized','outerposition',[0 0 1 1])
fig.PaperPosition = [0 0 30 12];
n_fig_events = 12;
j = 1 %-- specifies the page of events to represent
for i = (1+(j-1)*n_fig_events):j*n_fig_events
    subplot(2,6,i-n_fig_events*(j-1));
    surf(events(:,:,i),'FaceAlpha',0.99);
    view([90, 90]);
    axis([0 11 0 11 -200 1050])
    title({['Event: ', num2str(i) '/' num2str(n_events)]; ...
        [num2str(electrons(i)/1000) ' keV']; });
    xlabel('X(px)');
    ylabel('Y(px)');
    zlabel('ADU');
end

%% Histogram of events energies

figure('pos',[10 10 900 600]);
set(gca,'fontsize',20);
histogram(electrons/1000, 'Normalization', 'probability', 'BinWidth', 0.8);
hold on;
[count,bins]  = histcounts(electrons/1000);
title(['Number of electrons / event'],  'fontsize',20);
grid on;
xlabel('Number of collected electrons per event (thousands)', 'fontsize',20);
ylabel('Relative Frequency', 'fontsize',20);
