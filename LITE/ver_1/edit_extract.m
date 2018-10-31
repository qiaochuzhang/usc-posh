function [raw_data,iter]=edit_extract(inputs)
%edit all relevant schematics and extract data

%{
clear all; clc; clf; close all;
% load('testdata.mat');
load('usermain.mat');
tic;
%}

%inputs (combining editor and extractor inputs)
% %{
folder      = inputs{1};
schem       = inputs{2};
comps       = inputs{3};
prefixes    = inputs{4};
vals        = inputs{5};
iter        = inputs{6};
probe       = inputs{7};
dt          = inputs{8};
interpolate = inputs{9};
exe_path    = inputs{10};
%}

%orders of magnitude (using reciprocal convention, which treats prefixes as units)
femto=1e15; pico=1e12; nano=1e9; micro=1e6; milli=1e3; centi=1e2; 
kilo=1e-3; Mega=1e-6; Giga=1e-9; Tera=1e-12; Peta=1e-15;

%=========================================================================%
%SCHEMATIC EDITOR CODE====================================================%
%=========================================================================%

%define list of words for generating a new line in new schematic
NLW={'SHEET','WIRE','FLAG','SYMBOL','WINDOW','SYMATTR','TEXT','IOPIN'};
nnl=size(NLW,2); %number of exceptions for newline writing

%read in and store various subcircuit schematics/netlists as string
%cell structures:
fn=[folder,'/',schem,'.asc'];  
fid=fopen(fn); %generate file IDs
temp=textscan(fid,'%s')'; %convert to space-delimited string cell
temp=temp{1}; %organize scanned subcircuit schematic subcells

%number of components array
nc=size(comps,2);

%edit top-level schematic (no need to edit subcircuits!)
    
%check to see whether corresponding component value is a scalar (change
%.param value) or an array (change .param to .step param)
idi=zeros(1,nc);
ind=zeros(1,nc);
for i=1:nc %loop through all components in top-level schematic
    
    %detect whether component exists at all in subcircuit
    idi(i)=sum(ismember(temp,comps{i}));
    
    if idi(i)~=0 %(parameter exists)
        
        %determine index location for current schematic component
        ind(i)=find(ismember(temp,comps{i}));
        
        %check size of value cell
        if size(vals{i},2)==1 %change fixed parameter value
            
            valstr{i}=[num2str(vals{i}),prefixes{i}]; %convert values to string format
            
            %need to check whether directive is fixed or step
            %if fixed, do nothing but replace value
            %if stepped, change directive from .step param to .param
            
            %need to first find directive associated with current component
            test1=prod(ismember(temp{ind(i)-1},'!.param'));
            test2=prod(ismember(temp{ind(i)-2},'!.step'));
            if test1==1 %simply change value
                temp{ind(i)+1}=valstr{i}; %change parameter value (string format)
            end
            if test2==1 %change from .step param to .param
                temp{ind(i)-2}='!.param';
                temp{ind(i)-1}=' ';
                temp{ind(i)+1}=valstr{i}; %change parameter value (string format)
                temp{ind(i)+2}=' ';
                temp{ind(i)+3}=' ';
            end
            if test1~=1 && test2~=1
             %component is not associated with either .param nor .step param directive
                error([comps{i},' requires either a .PARAM or .STEP PARAM directive!']);                    
            end
            
        elseif size(vals{i},2)==3 %change .step param values
            
            tempSP=[];
            for k=1:3
                tempSP=[tempSP,' ',num2str(vals{i}(k)),prefixes{i}]; %convert values to string format
            end
            valstr{i}=tempSP;
                
            %need to check whether directive is fixed or step
            %if fixed, change directive from .param to .step param
            %if stepped, do nothing but replace [min,max,increment] values
            
            %need to first find directive associated with current component
            test1=prod(ismember(temp{ind(i)-1},'!.param'));
            test2=prod(ismember(temp{ind(i)-2},'!.step'));
            
            if test1==1 %change from .param to .step param
                temp{ind(i)-1}='!.step param';
                temp{ind(i)+1}=valstr{i}; %replace single value with [min,max,increment] string array
              end
            if test2==1 %change .step param values
                %change parameter value (string format)
                for k=1:3
                    temp{ind(i)+k}=[num2str(vals{i}(k)),prefixes{i}];
                end
            end
            if test1~=1 && test2~=1
                %component is not associated with either .param nor .step param directive
                error([comps{i},' requires either a .PARAM or .STEP PARAM directive!']);                    
            end            
        else %parameters must either be fixed (.param) or stepped (.step param)
            error(['"',comps{i},'" must either be a single ',...
                '.PARAM value or a .STEP PARAM array [min,max,increment].']);
        end        
    else        
        error(['Parameter "',comps{i},'" does not exist in schematic "',schem,'".']);        
    end    
end

%rewrite schematic file with updated value(s)
%only change top-level schematic
fid=fopen([folder,'/',schem,'_edit.asc'],'w');
fprintf(fid,[temp{1},' ',temp{2}]); %version #

for n=3:size(temp,1)
    for k=1:nnl %check 
        idx(k)=prod(ismember(temp{n},NLW{k})); %ID for New Line
    end
    if sum(idx)==1 %meaning that one of the key words was found
        fprintf(fid,'\n');
    end
    %         if i==libstr; for j=1:numel(test2{i}); if test2{i}(j)=='\'; test2{i}(j)='/'; end; end; end;
    fprintf(fid,[temp{n},' ']);
end
fprintf(fid,'\n');
fclose(fid);

iter = iter+1;

%=========================================================================%
%EXTRACTION CODE (UNGLAUB)================================================%
%=========================================================================%

%extract all relevant data from LTspice simulation
% File execution =========================================================%
exe_path=strtrim(exe_path);

%extract data from EDITED schematic (not original)
schem = [schem,'_edit'];
schemfile=[folder,'/',schem,'.asc']; %schematic filename
% netfile=[folder,'/',schem,'.net']; %corresponding netlist filename
rawfile=[folder,'/',schem,'.raw']; %raw filename generated by LTspice

%execution commands
% run_schem=['"' exe_path '" -ascii -b -run "' netfile '"'];
run_schem=['"' exe_path '" -ascii -b -run "' schemfile '"'];
system(run_schem);

%look for LTspice executable file
if ~exist(exe_path,'file')
    error('LTspice executable was not found. Either place it in the path or use the asc-based code.');
end

% Data extraction =======================================================%
%determine type of simulation
fid=fopen(rawfile);
temp1=textscan(fid,'%s')';
temp1=temp1{1};
temp2=find(ismember(temp1,'Variables:'));

novars=str2num(temp1{temp2(1)+1}); %number of variables (needed to extract values) 
nopts=str2num(temp1{temp2(1)+4}); %number of sample points
atype=temp1{temp2(2)+3}; %analysis type
nodeno=str2num(temp1{find(ismember(temp1,probe))-1}); %relevant node number (probe)
data_start=find(ismember(temp1,'Values:'));

if ismember(atype,'time') %.tran analysis
    
    %begin collecting time and probe values
    time=zeros(1,nopts); %preallocate time array
    vout=zeros(1,nopts); %preallocate output voltage array
    for i=1:nopts
        time(i)=str2num(temp1{data_start+2+(i-1)*(novars+1)});
        vout(i)=str2num(temp1{data_start+2+nodeno+(i-1)*(novars+1)});
    end    
    
    %time series interpolation (since LTspice doesn't use fixed timesteps)
    if interpolate==0
        time_interp = time; %new time series
        vout_interp = vout; %new voltage series
    elseif interpolate==1
        time_interp = time(1):dt:time(end); %new time series
        vout_interp = interp1(time,vout,time_interp); %new voltage series
    end
  
    %compute Fourier transform
    tlen=length(time_interp); %number of samples in time array
    dt_nano=mean(diff(time_interp*nano)); %average sampling period (s)->(ns)
    fs=1/dt_nano; %sampling frequency (GHz)
    fN=fs/2; %Nyquist frequency (GHz)
    freq=linspace(-fN,fN,tlen); %frequency array (GHz)
    FFT=fftshift(fftn(ifftshift(vout_interp)))/sqrt(numel(vout_interp)); %take FFT and shift
    
    %extracting nominal frequency of signal
%     disp('Computing output signal frequency...');
    Fout=abs(FFT); %compute magnitude
    mfo=max(Fout); %find peak
    Fout=Fout/mfo; %normalized FFT magnitude
    hf=round(length(freq)/2)+2; %ignore DC offset; examine positive half
    fo=max(Fout(hf:end));
    ftemp=find(Fout==fo);
    f0=freq(max(ftemp))*(Mega/Giga); %choose positive frequency (GHz)->(MHz)
    
    %append to output data structure
    raw_data.time = time_interp; %.tran sweep
    raw_data.vout = vout_interp; % signal measurement (V)
    raw_data.freq = freq;        % FFT frequency array (GHz)
    raw_data.fft  = FFT;         % FFT of signal (complex)
    raw_data.f0   = f0;          % nominal signal frequency (MHz)
    
elseif ismember(atype,'voltage') %.dc analysis
    
    %begin collecting DC and probe values
    dc=zeros(1,nopts); %preallocate time array
    vout=zeros(1,nopts); %preallocate output voltage array
    for i=1:nopts
        dc(i)=str2num(temp1{data_start+2+(i-1)*(novars+1)});
        vout(i)=str2num(temp1{data_start+2+nodeno+(i-1)*(novars+1)});
    end
    
    %append to output data structure
    raw_data.dc=dc;     %.dc sweep
    raw_data.vout=vout; % measurement

elseif ismember(atype,'frequency') %.ac analysis
    
    %begin collecting frequency and probe values
    ac=zeros(1,nopts); %preallocate time array
    vout=zeros(1,nopts); %preallocate output voltage array
    for i=1:nopts
        temp=str2num(temp1{data_start+2+(i-1)*(novars+1)});
        ac(i)=temp(1);
        temp=str2num(temp1{data_start+2+nodeno+(i-1)*(novars+1)});
        vout(i)=temp(1)+1i*temp(2); %voltage is generally complex
    end
    
    %append to output data structure
    raw_data.ac=ac;     %.ac sweep
    raw_data.vout=vout; % measurement
    
end

fclose('all'); %close any/all instances of LTspice
% toc

end