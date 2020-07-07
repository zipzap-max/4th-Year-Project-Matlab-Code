%% Loading Images
clear all
close all
file=uigetfile('.mat','Please select a Matlab workspace file containing your phantom datasets of label1, label2 and difference MR30 images.');
load(file)
TL=input('What is the label duration TL?');
flowrate=input('What is the Flowrate?');
Imagetype = menu('Is your image a brain or phantom image?','Brain','Phantom');
slicemin=3;
slicemax=7;
N=length(diff(:,1,1,1)); %number of pixels


%% Establishing what type of images to use

%Choosing a threshold level
if flowrate >= 200
    threshold =1000;
else
    threshold =500;
end

choice = menu('Do you want to use raw Label-Control image or Scanner processed MR30 image?','Label-Control','MR30');
if choice ==1
    Data=label1-label2;
    imagetype = 'Label-Control';
end
if choice ==2
    Data=diff;
    imagetype = 'MR30';
end

imagesc(Data(:,:,4,7));
colormap(hot)
colorbar
pbaspect([1 1 1 ])
%% Observing how the signal changes across all PLD's in slice 5
h= diff(:,:,5,:);
figure; montage(h,'DisplayRange',[0,1000]) %all PLDs
colormap(hot)
title(['Signal from slice 5'])
axis off
pbaspect([1 1 1])
colorbar
%% Applying estimated values for f, dt, Tau and fitting

%setting up the variables for the model
if Imagetype == 1
    lamda= 0.9;
else
    lamda = 0.32;
end

T1= 1.9; %seconds
T1b= 1.9; %seconds
dt= 0.5; %seconds
f= 60./6000; %perfusion
tau= 1;  %seconds
M0= 1;
truef=f;  
T1app =1./(1./T1+f./lamda);

%Real Data pCASL sequence
t= [0.05, 0.25, 0.5, 0.750, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5]; %seconds
x=1;
sliceno=1;
slicenomax=slicemax-slicemin+1;
Coordinates= zeros(slicenomax,2,(N).^2);
Norm = zeros(slicenomax, length(t), (N).^2);

for slice=slicemin:slicemax
    
    for n=1:N
    
        for m=1:N
        
    %thresholding useing 500 for low flow rates and 1000 for high flowrates
        if max(Data(n,m,slice,5:9))<=threshold
        else
            clean(1,:)=Data(n,m,slice,:);

            fort=clean./max(clean);
            for u = 1:length(t)
                a=fort(1,u);
                Norm(sliceno,u,x)=a;
            end
    %xcoordinates
            coordinates(x,1)=n;
    %ycoordinates
            coordinates(x,2)=m;
            x=x+1;
        end
        
        end
        
    end
    A = exist('coordinates');
    if A== 1
    Coordinates(sliceno,:,1:length(coordinates))=transpose(coordinates(:,:));
    sliceno=sliceno+1
    end
    clear coordinates
    x=1;
end
%% Flexfitting

dt= [0.2,1.5]; %seconds
f= [0.008,0.03]; %perfusion
tau= [0.9,1.1];  %seconds
t1= [0.05, 0.25, 0.5, 0.750, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5]; %seconds
t2= t1+ 1;

%Choosing which starting values to use
     T1app =1./(1./T1+f(1)./lamda);     
     fnew=f(1);
     dtnew=dt(1);
     taunew=tau(1);
     
clear sliceno

for sliceno = 1:slicenomax
    activevoxels = 0;
    while Coordinates(sliceno,1,activevoxels+1) ~= 0
        activevoxels=activevoxels+1;
    end
for n= 1:activevoxels;    
             clear y estimates
                     
                     y=Norm(sliceno,1:11,n);
                     sum = mean(y);
                     
if sum ~= 0                    
count=1;                     

                   
                     [estimates, model] = perfusionfit_3params(t1,y,M0,fnew,T1b,T1app,dtnew,taunew);

                     f1=estimates(1,:);
                     dt1=estimates(2,:);
                     tau1=estimates(3,:);            
%     Shifting the axis for negative dt
%     if f1==0
%         t2= t1+ TL;
%         [estimates, model] = perfusionfit_3params(t2,y,M0,fnew,T1b,T1app,dtnew,taunew);
%         f1=estimates(1,:);
%         dt1=estimates(2,:)-TL;
%         tau1=estimates(3,:);
%     end
                    fs(1,count)=f1;
                    dts(1,count)=dt1;
                    taus(1,count)=tau1;
    

                    count=count+1;
            end
    %takes the mean fit for each voxel for many values of perfusion, dt,
    %tau
    farray(sliceno,n)=mean(fs);
    dtarray(sliceno,n)=mean(dts);
    tauarray(sliceno,n)=mean(taus);
% end 
end
end
%% f mapping

N=length(diff(:,1,1,1));    
empty=zeros(N,N,slicenomax);
for sliceno = 1:slicenomax
for n=1:length(Coordinates(sliceno,1,:))
    if Coordinates(sliceno,:,n)== 0 
    else
    empty(Coordinates(sliceno,1,n),Coordinates(sliceno,2,n),sliceno)=farray(sliceno,n).*6000;
    end
end
end

figure; montage(empty,'DisplayRange',[0,500]) %all slices
colormap(hot)
title(['Perfusion Weighted Image of z= ', num2str(slicemin), ' to ', num2str(slicemax), ' slice through the capillary bed.'])
axis off
pbaspect([1 1 1])
colorbar
%% f histogram
clear y
y= (transpose((real(empty(:)))));
nbins= 100;
figure; histogram(nonzeros(y),nbins); 
title(['CBF Histogram flowrate ', num2str(flowrate), ' label duration ', num2str(TL)]);
hold on
graph1 = plot(fittedmodel);
set(graph1,'LineWidth',2);
legend('Histogram','Two Gaussian fit')
xlabel('CBF (ml/100g/min)')
ylabel('Voxels')
axis([-inf inf 0 inf])

mean_CBF = mean(nonzeros(abs(empty(:))))
stdev_CBF = std(nonzeros(abs(empty(:))))

mode_CBF = mode(nonzeros(empty(:)));
median_CBF = median(nonzeros(empty(:)));

clear cfx cfy
[N,edges]= histcounts(nonzeros(y),nbins);
binwidth= abs(edges(1)-edges(2));
for n=1:(length(N))
centres(n)= edges(n)+ binwidth/2;
end

% Clipping the histogram to use range -200:600
k=1;
for n=1:length(N)
    if centres(n) <350
        if centres(n) >-200
        cfy(k)= N(n);
        cfx(k)=centres(n);
        k=k+1;
        end
    end
end
%% dt mapping

N=length(diff(:,1,1,1));    
empty=zeros(N,N,slicenomax);
for sliceno = 1:slicenomax
for n=1:length(Coordinates(sliceno,1,:))
    if Coordinates(sliceno,:,n)== 0 
    else
    empty(Coordinates(sliceno,1,n),Coordinates(sliceno,2,n),sliceno)=dtarray(sliceno,n);
    end
end
end

figure; montage(empty,'DisplayRange',[0,1]) %all slices
colormap(hot)
title(['Delay Time Weighted Image of z= ', num2str(slicemin), ' to ', num2str(slicemax), ' slice through the capillary bed.'])
axis off
pbaspect([1 1 1])
colorbar
%% dt histogram
clear y
y= (transpose((real(empty(:)))));
nbins= 100;
figure; histogram(nonzeros(y),nbins); 
title(['dt Histogram flowrate ', num2str(flowrate), ' label duration ', num2str(TL)]);
hold on
% graph1 = plot(fittedmodel);
% set(graph1,'LineWidth',2);
legend('Histogram','Two Gaussian fit')
xlabel('Time (s)')
ylabel('Voxels')
axis([-inf inf 0 inf])



%% tau fitting

N=length(diff(:,1,1,1));    
empty=zeros(N,N,slicenomax);
for sliceno = 1:slicenomax
for n=1:length(Coordinates(sliceno,1,:))
    if Coordinates(sliceno,:,n)== 0 
    else
    empty(Coordinates(sliceno,1,n),Coordinates(sliceno,2,n),sliceno)=tauarray(sliceno,n);
    end
end
end

figure; montage(empty,'DisplayRange',[0,3]) %all slices
colormap(hot)
title(['Tau Weighted Image of z= ', num2str(slicemin), ' to ', num2str(slicemax), ' slice through the capillary bed.'])
axis off
pbaspect([1 1 1])
colorbar
%% Tau histogram
clear y
y= (transpose((real(empty(:)))));
nbins= 100;
figure; histogram(nonzeros(y),nbins); 
title(['Tau Histogram flowrate ', num2str(flowrate), ' label duration ', num2str(TL)]);
hold on
% graph1 = plot(fittedmodel);
% set(graph1,'LineWidth',2);
legend('Histogram','Two Gaussian fit')
xlabel('Time (s)')
ylabel('Voxels')
axis([-inf inf 0 inf])