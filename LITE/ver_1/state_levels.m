function [levels]=state_levels(signal)
%compute the state levels by taking averages...

%{
%testing the function...
npts=1000;
time=linspace(0,6*pi,npts);
duty=0.5;
noise=0.05;
signal=square(time)+noise*(rand(1,npts)-rand(1,npts));

figure(1)
plot(signal);

%generate signal histogram
hsig=histogram(signal,1000);

data=hsig.Data;
values=hsig.Values;
%}

%first find midpoint of signal levels
midpt=mean(signal);
%find average of lower levels
indlow=find(signal<midpt);
siglow=signal(indlow);
low=mean(siglow);
%find average of upper levels
indhigh=find(signal>=midpt);
sighigh=signal(indhigh);
high=mean(sighigh);

%export results
levels=[low,high];

end