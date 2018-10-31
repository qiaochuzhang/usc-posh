% User input script
clear all; clc; clf; close all;

%orders of magnitude (using reciprocal convention, which treats prefixes as units)
femto=1e15; pico=1e12; nano=1e9; micro=1e6; milli=1e3; centi=1e2; 
kilo=1e-3; Mega=1e-6; Giga=1e-9; Tera=1e-12; Peta=1e-15;

debug_mode = 1; %'1'=use preset values to run trial example

% Input target metrics -------------------------%

disp(['|=====================================================================|']);
disp(['|============================VCO OPTIMIZER============================|']);
disp(['|=====================================================================|']);

disp(' ');

debug_choice=input(['Choose to run optimizer with preset parameters (1) or manually select (2): ']);

if debug_choice==1
    
    NS = 9; %number of inverter stages in VCO circuit
    FT = 337; %target frequency (MHz)
    
    nvars = 3; %number of variable parameters
    nsteps = 6; %number of steps along a variable parameter axis for seed grid generation
    vbnds = [500,1000;500,1000;5,15]; %[min,max] variable bounds
    dvar = [10,10,0.065]; %step size for variables (used in gradient descent)
    
    %evaluation bounds for validity of simulation results
    SLB = [0.2,0.9]; %[min,max] state level bounds (V)
    DCB = [0.35,0.65]; %[min,max] duty cycle bounds
    %search tolerances
    eps = 10; %percent distance of target frequency to search in database (%)
    tol = 0.5; %percent tolerance for gradient/steepest descent (%)
    nit = 3; %number of iterations to scan for convergence to incorrect minimum
    
    %DEFAULT CIRCUIT PARAMETER VALUES    
    
    %list of all top-level parameters
    Tmax = 400; %total simulation time for VCO circuit (ns)
    dt   = 10; %maximum simulation timestep/interpolation timestep (ps)
    Vs   = 1.1; %source voltage (V)
    Vctr = 0.29; %control voltage (V)
    %list of all fixed subcircuit parameters
    Lp = 65; %pmos transistor length (nm)
    Ln = 65; %pmos transistor width (nm)
    R = 10; %bias resistor (kOhm)    
    %user-chosen (arbitrary) variable values (avoid .step param for now)
    Wp  = 300; %(nm) 
    Wn  = 200; %(nm)
    Wnb = 1; %(um)
    
    %variable parameter bounds & increment sizes (usage: [min,max,increment])
    % Wp  = [300,1800,1500]; %(nm) 
    % Wn  = [200,1200,1000]; %(nm)
    % Wnb = [1,21,20]; %(um)
    
    % Wp  = [300,1800,150]; %(nm) 
    % Wn  = [200,1200,100]; %(nm)
    % Wnb = [1,21,2]; %(um)  
    
elseif debug_choice==2
    
    disp(['USER INPUTS']);
    disp([' ']);
    disp(['Circuit: VCO']);
    input(['Choose technology file (only "PTM" available for now): '],'s');
    input(['Pick model size (nm) (only "65" available for now): ']);
    NS=input(['Pick number of inverter stages, N (only "9" available for now): ']);
    
    disp([' ']);
    disp(['Target metrics']);
    FT=input(['Pick target nominal frequency (MHz), f0 = ']);
    
    disp([' ']);
    disp(['Circuit parameter values']);
    PPD=input(['Use default values for all non-variable circuit parameters? (1=Yes,0=No): ']);
    
    
    %Parameter Prompt Decision
    if PPD==1    
        nvars=3;
        vbnds = [150,6*150;100,6*100;1,6*1]; %[min.max] variable bounds
        dvar = [150,100,1]*0.1; %step size for variables (used in gradient descent)
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
    elseif PPD==0        
        nvars=input(['Choose number of variable parameters: ']);
        for i=1:nvars
            vbnds(i,:)=input(['Choose [min,max] bounds for variable ',num2str(i),': ']);
            dvar(i)=input(['Choose gradient descent step size for variable ',num2str(i),': ']);
        end
        %list of all top-level parameters
        Tmax = input(['Choose total cicuit simulation time Tmax (ns): ']);
        dt   = input(['Choose interpolation timestep dt (ps): ']);
        Vs   = input(['Choose source voltage Vs (V): ']);
        Vctr = input(['Choose control voltage Vctr (V): ']);
        %list of all fixed subcircuit parameters
        Lp   = input(['Choose PMOS transistor length (nm): ']);
        Ln   = input(['Choose NMOS transistor length (nm): ']);
        R    = input(['Choose value of bias resistor R (kOhm): ']);
        %initial (arbitrary) variable values
        Wp   = input(['Choose a starting value for 1st variable, Wp (nm): ']);
        Wn   = input(['Choose a starting value for 2nd variable, Wn (nm): ']);
        Wnb  = input(['Choose a starting value for 3rd variable, Wn,b (um): ']);      
    end
    
    disp([' ']);
    disp(['Constraints']);
    SLB=input(['Choose [min,max] state level bounds for noise margin: ']);
    DCB=input(['Choose [min,max] duty cycle bounds: ']);
    eps=input(['Choose percent distance of target frequency for database search (%): ']);
    tol=input(['Choose percent tolerance for gradient/steepest descent (%): ']);
    nit=input(['Choose number of iterations to scan for convergence to incorrect minimum: ']);
    
end

%SWITCHES=================================================================%
plotting = 1; %option to generate plots throughout and/or at end of pipeline (default=0)
savefigs = 1; %option to automatically save generated figures
savedata = 1; %option to save all data after complete pipeline execution (default=0)
%switch array for database seed generator
switches=[plotting,savefigs,savedata];
interpolate = 1; %option to interpolate time series extracted from LTspice (default=1)

%METRICS==================================================================%
% Save user data
metrics{1} = FT;
metrics{2} = SLB;
metrics{3} = DCB;
metrics{4} = [eps,tol,nit];
metrics{5} = nvars;
metrics{6} = vbnds;
metrics{7} = dvar;
save(['metrics.mat'],'metrics');

%INPUTS===================================================================%

%executable paths for LTspice
lts_path=['C:\Program Files\LTC\LTspiceXVII\'];
exe_path=[lts_path,'XVIIx64.exe'];

% SCHEMATICS (TOP LEVEL BLOCK AND SUBCIRCUITS)
folder = 'VCO_files'; %location of schematic files
probe = 'V(osc)'; %voltage node data to extract
schem = ['VCO',num2str(NS)]; %top level block schematic
iter = 0; %initialize iteration index
%preparation of circuit parameter cell arrays (labels and values)
comps = {'Vs','Vctr','Tmax','Lp','Ln','R','Wp','Wn','Wnb'};
prefixes = {'','','n','n','n','k','n','n','u'};

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

%=========================================================================%

if savedata==1
    save(['usermain.mat']);
end

% %{

%must check whether database file exists or not
DB_exist_check = exist('database.mat','file');
MG_exist_check = exist('metric_grid.mat','file');
VA_exist_check = exist('var_arrays.mat','file');

if DB_exist_check~=2 %if no database is found in current folder
    disp(['Database not found! Look for metric grid for database seeding...']);
    if MG_exist_check==2 && VA_exist_check==2 %no database found but a metric grid found
        disp(['Metric grid found, proceeding to generate an initial database...']);
        %map valid grid into a rudimentary matrix database
        metric_grid=importdata('metric_grid.mat');
        var_arrays=importdata('var_arrays.mat');
        [database]=database_gen(metric_grid,var_arrays);
    else %no database OR metric grid found
        disp(['No metric grid found! Proceeding to construct metric grid for database seeding...']);
        nsteps=input(['Choose number of steps along each variable axis for database initialization N = ']);
        metrics{8} = nsteps;        
        %generate grid of valid metrics
        [metric_grid,var_arrays]=grid_gen(switches,inputs,metrics,nsteps);
        %map valid grid into a rudimentary matrix database
        [database]=database_gen(metric_grid,var_arrays);
    end
end

% Other metrics that could be additionally considered:
% -VCO tuning range (bandwidth)
% -Spectral purity (timing jitter, phase noise)
% -Power consumed
% -Load pulling
% -VCO sensitivity

% Priority weighting array for different metrics will also be incorporated 
% into a future version.

% Initiate optimization pipeline
pipeline(switches,metrics,inputs);

%load all generated data files
% load('var_arrays.mat');
% load('metric_grid.mat');
% load('seeddata.mat');
load('database.mat');
load('metrics.mat');
load('pipeline_data.mat');

%close any/all LTspice instances running in background
fclose('all');


