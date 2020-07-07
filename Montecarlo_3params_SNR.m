clear, close all
%% Establishing starting values
ff= [350,300,200,100];
dt= 0.2; %seconds
lamda= 0.34;
T1= 1.8; %seconds
T1b= 2.1; %seconds
M0= 1;
tau= 1;  %seconds

repeats = input('How many repeat images were taken?');
percentgeerror=1./sqrt(repeats - 1);
%% Running a loop for multiple flowrates u

for u=1:4
t=linspace(0,3); %seconds
f= ff(u)./6000; %perfusion

%pCASL Sequence

T1app =1./(1./T1+f./lamda);

T1r =1./(1./T1app-1./T1b);

%Simulating magnetization at flowrate ff(u)
for n=1:length(t)
    
if t(1,n) < dt
        M(1,n)= 0;
else
    if t(1,n) < dt + tau
        M(1,n)= 2.*100.*M0.*f.*T1app.* exp(-dt./T1b).*(1-exp(-(t(1,n)-dt)./T1app));
    else
        M(1,n)= 2.*100.*M0.*f.*T1app.* exp(-dt./T1b).*exp(-(t(1,n)-dt-tau)./T1app).*(1-exp(-tau./T1app));
    end
end

end

% plot(t,M)
% title('pCASL Sequence')
%% 
%Adding random noise to the data points then fitting for 1000 loops

hold on 
t= [0.05, 0.25, 0.5, 0.750, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5];
peak= max(M);
clear M
%MonteCarlo
for m=1:1000 %% Number of repeats in montecarlo = 1000

for n=1:length(t)
    
if t(1,n) < dt
        Mag(u,n)= 0;
else
    if t(1,n) < dt + tau
        Mag(u,n)= 2.*100.*M0.*f.*T1app.* exp(-dt./T1b).*(1-exp(-(t(1,n)-dt)./T1app)); %%This is as a percentage change
    else
        Mag(u,n)= 2.*100.*M0.*f.*T1app.* exp(-dt./T1b).*exp(-(t(1,n)-dt-tau)./T1app).*(1-exp(-tau./T1app)); %%This is as a percentage change
    end
end
%

error=percentgeerror*(-1+2*rand);
M(1,n)= Mag(u,n)+error;


end

%fminsearch fitting


    [estimates, model] = perfusionfit_3params(t,M,M0,f,T1b,T1app,dt,tau);

   perf(1,m)=estimates(1,:);
   deltat(1,m)=estimates(2,:);
   tau_time(1,m)=estimates(3,:);
   
    

end


% perf(1,m)=estimates(1,:);
% deltat(1,m)=estimates(2,:);
% tau_time(1,m)=estimates(3,:);

perfusion= mean(perf);
perfusionunits=perfusion*6000;
error_f= std(perf);
error_funits=error_f*6000;
percentageerror(u,:) = (error_f./perfusion);
dt=mean(deltat);
error_dt=std(deltat);
tau=mean(tau_time);
error_tau=std(tau_time);
disp(['Perfusion is measured as ',num2str(perfusionunits),' ml/100g/min with an error of ',num2str(error_funits),' ml/100g/min.'])
disp(['dt is measured as ',num2str(dt),' s with an error of ',num2str(error_dt),' s.'])
disp(['Tau is measured as ',num2str(tau),' s with an error of ',num2str(error_tau),' s.'])
t=linspace(0,3);

for n=1:length(t)
    
if t(1,n) < dt
        M(1,n)= 0;
else
    if t(1,n) < dt + tau
        M(1,n)= 2.*100.*M0.*perfusion.*T1app.* exp(-dt./T1b).*(1-exp(-(t(1,n)-dt)./T1app));
    else
        M(1,n)= 2.*100.*M0.*perfusion.*T1app.* exp(-dt./T1b).*exp(-(t(1,n)-dt-tau)./T1app).*(1-exp(-tau./T1app));
    end
end

end
Marray(u,:)=M;
end
%% Plotting the simulated perfusion curves with error due to monte carlo STD
t1= [0.05, 0.25, 0.5, 0.750, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5];
Marraynew= Marray./max(max(Marray));
%%Magg= Mag./max(max(Marray));
for n=1:4
plot(t,Marray(n,:))
hold on
error1= Mag(n,:).*percentageerror(n);
errorbar(t1,Mag(n,:),error1,'*')
end
legend('350 flowrate','','300 flowrate','','200 flowrate','','100flowrate','')
title('Monte Carlo pCASL Sequence with a percentage error of ', percentageerror)
xlabel('Time (s)')
ylabel('{\Delta} M')

