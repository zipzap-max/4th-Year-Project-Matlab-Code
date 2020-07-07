close all; clear all; clc;
%% load the label, control and scanner generated difference images

file=uigetfile('.mat','Please select a Matlab workspace file containing your phantom datasets of label1, label2 and difference MR30 images.');
load(file)
dim = size(label1);
threshold = 500;
Imagetype = menu('Is your image a brain or phantom image?','Brain','Phantom');


%% show image of label 2 to check image quality

h= diff(:,:,:,8);
figure; montage(h,'DisplayRange',[0,1000]) %all slices
colormap(hot)
title(['Signal from slice 3-7'])
axis off
pbaspect([1 1 1])
colorbar
figure; montage(label2(:,:,:,11),'DisplayRange',[0,50000]) %all slices
title('Control image, PLD = 2500ms')
colorbar


%% only use central 5 slices
%Checking at which central five slices the porous material is in
%(currently 3-7)
slicemin=3;
slicemax=7;

label1 = label1(:,:,slicemin:slicemax,:);
label2 = label2(:,:,slicemin:slicemax,:);
diff = diff(:,:,slicemin:slicemax,:);
M0_image = label2(:,:,:,8);

%% mask option 1 -  hough transform (ish)

% First using matlab functions to define circle
% Only performing Hough Transform if using phantom data

if Imagetype ==2;
image = label2(:,:,3,8);
[centers, radii, metric] = imfindcircles(image,[20 500]);
radii = radii+2;
% h = viscircles(centers(1,:),radii);
% figure; imagesc(label2(:,:,3,8)); hold on; viscircles(centers(1,:),radii);

[X,Y] = meshgrid(1:size(image,1), 1:size(image,2));
circle_mask = sqrt(((X-centers(1)).^2)+((Y-centers(2)).^2)) <= radii;
% figure; imagesc(circle_mask);
% title('Circle mask')

% apply mask to label 2 TI 2500 to get M0


M0_circle_masked=M0_image.*circle_mask;
median_circle = median(nonzeros(M0_circle_masked));
end

%% mask option 2 - thresholding
image = label2(:,:,:,8);
M0_image = label2(:,:,:,8);
dim_mask = size(image);
thresh_mask = ones(dim_mask);
for i = 1:dim_mask(1)
    for j = 1:dim_mask(2)
        for k = 1:dim_mask(3)
            
            if image(i,j,k) <threshold
                thresh_mask(i,j,k) = 0;
            else
            end
        end
    end
end

% figure; montage(thresh_mask,'DisplayRange',[0,1]) %all slices
M0_thresh_masked=M0_image.*thresh_mask;
median_thresh = median(nonzeros(M0_thresh_masked));

%% mask option 3 - from scanner diffs

maskoff = diff(:,:,:,8);
dim = size(label2);
for p = 1:dim(3)
    for n=1:dim(2)
        for m=1:dim(1)
            if maskoff(n,m,p) < 200
                mask2(n,m,p) = 0;
            else
                mask2(n,m,p) = 1;
            end
        end
    end
end

% figure; montage(mask2,'DisplayRange',[0,1]) %all slices
M0_mask2 = M0_image.*mask2;
median_mask2 = median(nonzeros(M0_mask2));

%% choosing which mask to use

choice = menu('Which mask do you want to use?','Hough','Threshold','MR30');
if choice == 1
    M0 = median_circle;
    mask = circle_mask;
    masktype= 'Hough';
else
    if choice == 2
        M0 = median_thresh;
        mask = thresh_mask;
        masktype=  'Threshold';
    else
        M0 = median_mask2;
        mask = mask2;
        masktype= 'MR30';
    end
end

%% Applying estimated values for f, dt, Tau and fitting

%setting up the variables for the model
%If Brain lamda = 0.9, if phantom lamda = 0.32
if Imagetype == 1
    lamda= 0.9;
else
    lamda = 0.32;
end
T1= 1.9; %seconds
T1b= 1.8; %seconds
dt= 0.9; %seconds
f= 60./6000; %perfusion
tau= 1.8;  %seconds
alpha = 0.85; %pCASL labelling efficiency
PLD = 1.750;
scan = 'pCASL';

%     scan= input('Is scan pCASL or PASL?','s');
%     while scan ~= 'pCASL' or 'PASL'
%         scan= input('Please enter a valid scan type pCASL or PASL?','s');
%     end


for p=1:dim(3)
    for n= 1:dim(2)
        for m = 1:dim(1)
            
            if scan == 'pCASL'
                CBF(n,m,p) = (6000.*lamda.*(label1(n,m,p,8) - label2(n,m,p,8)).* exp(PLD./T1b))./ (2.* alpha.*T1b.*M0.*(1 - exp(-tau./T1b))); %%check this should be last PLD
            end
            
            %     if scan == 'PASL'
            %         CBF(n,m) = (6000.*lamda.*(label1(n,m,5,8) - label2(n,m,5,8)).* exp(PLD./T1b))./ (2.* alpha.*T1b.*label2(n,m,5,11));
            %     end
            
             if  CBF(n,m,p) > 1000
                 CBF(n,m,p) = 0;
             end
            
            if isnan(CBF(n,m,p));
                CBF(n,m,p) = 0;
            end
        end
    end
end
%% Applying the CBF mask

%figure; imagesc(CBF(:,:,5)) %example slice
figure; montage(CBF,'DisplayRange',[0,500]) %all slices
colormap(hot)
title('Simple CBF mapping Flowrate 350 TL 1800')
axis off
pbaspect([1 1 1])
colorbar%% mask CBF image

CBF_masked = CBF.*mask;
%figure; imagesc(CBF(:,:,4)) %example slice
title('Example CBF image')
figure; montage(CBF_masked,'DisplayRange',[0,500]) %all slices
colormap(hot)
colorbar
pbaspect([1 1 1])
title('CBF masked')
y= (transpose(real(CBF_masked(:))));

%% Plotting a histogram of the CBF distribution

nbins= 100;
figure; histogram(nonzeros(y),nbins); 
title(['CBF Histogram ', masktype]);
axis([-400 1000 -inf inf])
hold on
% graph1 = plot(fittedmodel);
% set(graph1,'LineWidth',3);
axis([-400 1000 -inf inf])
legend('Histogram','Two Gaussian fit')
xlabel('CBF (ml/100g/min)')
ylabel('Voxels')
%% mean etc...

mean_CBF = mean(nonzeros(abs(CBF_masked(:))))
stdev_CBF = std(nonzeros(abs(CBF_masked(:))))

mode_CBF = mode(nonzeros(CBF_masked(:)));
median_CBF = median(nonzeros(CBF_masked(:)));

%% Converting the histogram into data points 
%Can use cfx and cfy in cftool

clear cfx cfy
[N,edges]= histcounts(nonzeros(y),nbins);
binwidth= abs(edges(1)-edges(2));
for n=1:(length(N))
centres(n)= edges(n)+ binwidth/2;
end

% Clipping the histogram to use range -200:600
k=1;
for n=1:length(N)
    if centres(n) <1000
        if centres(n) >-600
        cfy(k)= N(n);
        cfx(k)=centres(n);
        k=k+1;
        end
    end
end