% function [metric_grid,var_arrays]=grid_gen(switches,inputs,metrics,nsteps)
%grid generator (precedes database seed generator)

%=========================================================================%
% %{
clear all;clc;clf;close all;

%orders of magnitude (using reciprocal convention, which treats prefixes as units)
femto=1e15; pico=1e12; nano=1e9; micro=1e6; milli=1e3; centi=1e2; 
kilo=1e-3; Mega=1e-6; Giga=1e-9; Tera=1e-12; Peta=1e-15;

plotting = 1; %option to generate plots throughout and/or at end of pipeline (default=0)
savefigs = 1; %option to automatically save generated figures
savedata = 1; %option to save all data after complete pipeline execution (default=0)
%switch array for database seed generator
switches=[plotting,savefigs,savedata];

NS = 9; %number of inverter stages in VCO circuit
FT = 345; %target frequency (MHz)

%bounds and number of steps (crucial parameters)
nvars = 3; %number of variable parameters
nsteps = 6; %number of steps along a variable parameter axis for seed grid generation
vbnds = [500,1000;500,1000;5,15]; %[min,max] variable bounds
dvar = [150,100,1]*0.1; %step size for variables (used in gradient descent)

%executable paths for LTspice
lts_path=['C:\Program Files\LTC\LTspiceXVII\'];
exe_path=[lts_path,'XVIIx64.exe'];
interpolate = 1; %option to interpolate time series extracted from LTspice (default=1)
folder = 'VCO_test'; %location of schematic files
schem = ['VCO',num2str(NS)]; %top level block schematic
%preparation of circuit parameter cell arrays (labels and values)
comps = {'Vs','Vctr','Tmax','Lp','Ln','R','Wp','Wn','Wnb'};
prefixes = {'','','n','n','n','k','n','n','u'};
probe = 'V(osc)'; %voltage node data to extract
iter = 0; %initialize iteration index

%evaluation bounds for validity of simulation results
SLB = [0.2,0.9]; %[min,max] state level bounds (V)
DCB = [0.35,0.65]; %[min,max] duty cycle bounds
%search tolerances
eps = 5; %percent distance of target frequency to search in database (%)
tol = 1; %percent tolerance for gradient/steepest descent (%)
nit = 5; %number of iterations to scan for convergence to incorrect minimum

%list of all top-level parameters
Tmax = 200; %total simulation time for VCO circuit (ns)
dt   = 10; %maximum simulation timestep/interpolation timestep (ps)
Vs   = 1.1; %source voltage (V)
Vctr = 0.29; %control voltage (V)        
%list of all fixed subcircuit parameters
Lp   = 65; %pmos transistor length (nm)
Ln   = 65; %pmos transistor width (nm)
R    = 10; %bias resistor (kOhm)   
%initial (arbitrary) variable values
Wp   = 300; %(nm) 
Wn   = 200; %(nm)
Wnb  = 1; %(um)   

%top-level
vals{1} = Vs;
vals{2} = Vctr;
vals{3} = Tmax;
%subcircuit constants
vals{4} = Lp;
vals{5} = Ln;
vals{6} = R;
%subcircuit variables
vals{7} = Wp;
vals{8} = Wn;
vals{9} = Wnb;

%define input cell array for schematic edit/data extraction
inputs{1}=folder;
inputs{2}=schem;
inputs{3}=comps;
inputs{4}=prefixes;
inputs{5}=vals;
inputs{6}=iter;
inputs{7}=probe;
inputs{8}=dt/pico; %MUST CONVERT TIME INCREMENT FROM (ps)->(s)!!!
inputs{9}=interpolate;
inputs{10}=exe_path;



% Save user data
metrics{1} = FT;
metrics{2} = SLB;
metrics{3} = DCB;
metrics{4} = [eps,tol,nit];
metrics{5} = nvars;
metrics{6} = vbnds;
metrics{7} = dvar;

%}
%=========================================================================%

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
    
%     [raw_data,~]=edit_extract1(inputs);
    [raw_data,~]=edit_extract2(inputs);
        
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
    save(['var_arrays.mat'],'var_arrays','-v7.3');
    save(['metric_grid.mat'],'metric_grid','-v7.3');
end

if plotting==1
    save('seeddata.mat');
    plotter_grid(switches,nsteps); %generate figures
    close all; %close all currently open figures
end


% end