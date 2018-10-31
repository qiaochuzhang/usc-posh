function [metric_grid,var_arrays]=grid_gen(switches,inputs,metrics,nsteps)
%grid generator (precedes database seed generator)

%{
clear all;clc;clf;close all;
load('testdata.mat');
plotting = 0; %option to generate plots throughout and/or at end of pipeline (default=0)
savefigs = 1; %option to automatically save generated figures
savedata = 1; %option to save all data after complete pipeline execution (default=0)
%switch array for database seed generator
nsteps = 5;
vbnds = [150,750;100,500;1,5];
tic;
%}

%switches array
plotting = switches(1);
savedata = switches(3);

%relevant metrics
SLB   = metrics{2}; %state level bounds (V)
DCB   = metrics{3}; %duty cycle bounds
nvars = metrics{5}; %number of variable parameters
vbnds = metrics{6}; %variable bounds

%variable parameter bounds & increment sizes (usage: [min,max,increment])
%generate arrays
for i=1:nvars
    var_arrays{i}=linspace(vbnds(i,1),vbnds(i,2),nsteps);
end

for k=1:nsteps
    disp(['k = ',num2str(100*(k-1)/nsteps),'%']);
for j=1:nsteps
    disp(['j = ',num2str(100*(j-1)/nsteps),'%']);
for i=1:nsteps
%     disp(['i = ',num2str(100*(i-1)/nsteps),'%']);
    
    %define single values for iteration loop
    Wp  = var_arrays{1}(i);
    Wn  = var_arrays{2}(j);
    Wnb = var_arrays{3}(k);

    %update input cell for editing/extraction
    vals{1} = inputs{5}{1};
    vals{2} = inputs{5}{2};
    vals{3} = inputs{5}{3};
    vals{4} = inputs{5}{4};
    vals{5} = inputs{5}{5};
    vals{6} = inputs{5}{6};    
    vals{7} = Wp;
    vals{8} = Wn;
    vals{9} = Wnb;
    inputs{5} = vals;
    
    %edit schematic and extract output signal
    [raw_data,~]=edit_extract(inputs);
        
    vout = raw_data.vout;
    f0   = raw_data.f0;
        
    % disp('Generating objective function grid...');
    
    data_in{1} = vout;
    reqs{1} = SLB;
    reqs{2} = DCB;
    
    %decide whether point is valid
    [td]=evaluator(data_in,reqs);
    
    if td==1
        metric_grid(i,j,k)=f0; %valid (units of MHz)
    elseif td==0
        metric_grid(i,j,k)=NaN; %NOT valid
    end
    
    fclose('all'); %close current LTS instance
    
end
end
end

%save out results
if savedata==1
    save(['var_arrays.mat'],'var_arrays');
    save(['metric_grid.mat'],'metric_grid');
end

if plotting==1
    save('seeddata.mat');
    plotter_grid(switches,nsteps); %generate figures
    close all; %close all currently open figures
end


end