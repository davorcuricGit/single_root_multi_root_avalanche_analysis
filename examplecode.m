clear all
load('ImgF_sample.mat')
%%

downsamplesize = 2;

load('sampleMask.mat');
maskOG = Mask; clear Mask
maskOG = maskOG/max(max(maskOG));

mask= spatialBlockDownsample(single(maskOG),downsamplesize, maskOG);
mask(mask~=0) = 1;
sz = size(maskOG) / downsamplesize;
validPixels = find(mask == 1);

threshold = 1;
ImgF2 = ImgF > threshold;

%% the radius was determined from the partial correlation function decay as a function of distance amongst pixels
radius = 16/downsamplesize;

%get the network
[~,network] = distance_network(sz(1),validPixels, radius);

%%
segments = get_non_zerosegments(ImgF2);

%% extract the avalanches for each non-zero segment, storing their initiation and termination points in all_clusters

%
%for each segment of the recording find clusters
SE = get_strel_dist(radius);
all_clusters = [];
next_id_offset = 0;
for s = 1:length(segments)
    % 1) dilate each frame

    seg = ImgF2(:,segments{s});
    seg = embeddIntoFOV(seg, validPixels, sz);

    f = imdilate(seg, strel('disk', radius));
    f = f.*mask;


    % 2) label clusters
    CC = bwconncomp(f, 6);
    labels3D = labelmatrix(CC);     % H×W×T of cluster IDs
    labeled = single(labels3D).*single(seg); %remove the dilation
    labeled = reshape(labeled, prod(sz), length(segments{s}));
    labeled = labeled(validPixels, :);

    clusters = getClusterStats(labeled, network);
    
    % Offset IDs to avoid label collisions as well as starttimes
    for i = 1:numel(clusters)
        clusters(i).id = clusters(i).id + next_id_offset;
        clusters(i).start_time = clusters(i).start_time + segments{s}(1);
    end

    all_clusters = [all_clusters, clusters];
    next_id_offset = max([all_clusters.id]);
end

clear s next_id_offset CC SE clusters labeled labels3D f
%%
S = [all_clusters.size];
D = [all_clusters.duration];
M = [all_clusters.merge_flag];

S_mergeless = S(M == 0);
D_mergeless = D(M == 0);

S_merged = S(M == 1);
D_merged = D(M == 1);

clf
X = {S, S_mergeless, S_merged};
exponent = [1.25, 1.5, 1.];
titles = {'All', 'Mergeless', 'Merged'};
for i = 1:3
    subplot(3,1,i)
[md, fq] = lnbin(X{i},30);
loglog(md, fq)
hold on
loglog(md, 0.5*md.^-exponent(i), '--k', 'DisplayName', ['\tau = ' num2str(exponent(i))])

title(titles{i})

if i == 2
ylabel('p(S)')
end

if i == 3
xlabel('Avalanche Size, S')
end
legend(legendUnq())
set(gca,'fontsize', 13)
end




%% Example: plotting root maps

density = zeros(sz);
for i =1:length(all_clusters)
    frame = zeros(1, length(validPixels));
    frame([all_clusters(i).roots]) = 1;
    frame = embeddIntoFOV(frame', validPixels, sz);
    density = density + frame;
end

imagesc(density)







function segments = get_non_zerosegments(ImgF)
%segment the recording into chunks of non-zero activity;
trace = sum(ImgF);
nanidx = find(isfinite(trace) == 0);
trace(nanidx) = 0;
[~,segments] = segmentSeries(trace);
segments = {segments{2:2:end}};

end




