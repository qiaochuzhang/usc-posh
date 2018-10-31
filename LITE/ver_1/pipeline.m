function pipeline(switches,metrics,inputs)
%main pipeline function (executes all auxiliary functions)
% clear all;clc;clf;close all;
TPT_tic=tic; %total pipeline time
appendDB=0; %switch to append KGD to database (default must be 0 to avoid copies)
disp('Taking in user specs and preparing circuit arrays...');

%{
load(['usermain.mat']);
%}

%extract target metrics and tolerances
FT    = metrics{1}; %target frequency (MHz)
SLB   = metrics{2}; %state level bounds (V)
DCB   = metrics{3}; %duty cycle bounds
eps   = metrics{4}(1); %percent distance of target frequency to search in database (%)
tol   = metrics{4}(2); %percent tolerance for gradient/steepest descent (%)
nit   = metrics{4}(3); %number of iterations to scan for convergence to incorrect minimum
nvars = metrics{5}; %number of variables
vbnds = metrics{6}; %variable parameter bounds
dvar  = metrics{7}; %variable step size for gradient descent

%define input cell array for schematic edit/data extraction
folder      = inputs{1};
schem       = inputs{2};
comps       = inputs{3};
prefixes    = inputs{4};
vals        = inputs{5};
iter        = inputs{6};
probe       = inputs{7};
dt          = inputs{8}; %already in units of (s)
interpolate = inputs{9};
exe_path    = inputs{10};

%print summary of user-input metrics/parameters to screen
disp([' ']);
disp(['Metrics:  Target nominal frequency f0 = ',num2str(FT),' MHz']);
disp(['          Noise margins: <=',num2str(SLB(1)),' and >=',num2str(SLB(2)),' V']);
disp(['          Duty cycle bounds: >=',num2str(DCB(1)),' and <=',num2str(DCB(2))]);
disp(['          Percent distance of target frequency for database search: ',num2str(eps),'%']);
disp(['          Local minimum tolerance for gradient descent search: ',num2str(tol),'%']);
disp(['          Number of Nesterov search iterations before giving up: ',num2str(nit)]);

savedata = switches(3);
%other switch values (only change for troubleshooting purposes)
pipe_switch = 1; %initiate data fitting & convexification (default=1)
%('1'=begin with database search,'2'=skip to gradient descent,'3'=skip to convex optimization)
localfail   = 0; %condition required to switch to data fitting subroutines (default=0)

%orders of magnitude (using reciprocal convention, which treats prefixes as units)
femto=1e15; pico=1e12; nano=1e9; micro=1e6; milli=1e3; centi=1e2; 
kilo=1e-3; Mega=1e-6; Giga=1e-9; Tera=1e-12; Peta=1e-15;

disp([' ']);
disp(['================================================================']);
disp(['INITIATING CIRCUIT OPTIMIZATION PIPELINE']);
disp(['================================================================']);
disp(['PHASE ONE: DATABASE SEARCH']);
P1T_tic=tic; %phase 1 timer start

if pipe_switch==1 %search database
    
    %[database search code goes here]
    [candidates] = database_search(metrics);
    NC = size(candidates,1); %(candidates are sorted from least to greatest frequency)
    
    if size(candidates)>0 %at least 1 candidate found
        disp(['Candidates found in database; test top candidate for tolerance criteria...']);
        
        %compute relative error
        relerr=100*abs(candidates(1,1)-FT)/FT;
                
        if relerr<=tol            
            disp([' ']);
            disp(['Top candidate falls within tolerance bounds; a solution has been found!']);
            disp(['Good design solution metric:']);
            final_f0=candidates(1,1);
            disp(['f0 = ',num2str(final_f0),' MHz']);
            final_err=relerr;
            disp(['Relative error = ',num2str(final_err),'%']);
            disp(['Good design solution parameters:']);
            %will want to automate this using the header in the future!!
            final_vars=[candidates(1,2),candidates(1,3),candidates(1,4)];
            disp(['Wp = ',num2str(final_vars(1)),' nm']);
            disp(['Wn = ',num2str(final_vars(2)),' nm']);
            disp(['Wn,b = ',num2str(final_vars(3)),' um']);            
            
            pipe_switch=0; %stop the flow and directly output circuit (do NOT append to database)
            solution=1; %a solution has been found; proceed to output the circuit
        else
            disp(['Top candidate outside tolerance bounds; begin local optimization phase...']);
            pipe_switch=2;
        end
        
        
    else %no candidates found
        disp(['No candidates found in database; perform data fitting and convex optimization.']);
        %no such solution exists in database; begin analysis pipeline
        pipe_switch=3;    
    end

end

P1T=toc(P1T_tic);
disp([' ']);
disp(['Phase 1 complete. Time = ',num2str(P1T),' seconds.']);

if pipe_switch==2 %gradient descent search
    disp(['----------------------------------------------------------------']);
    disp(['PHASE TWO: LOCAL SEARCH VIA GRADIENT DESCENT']);
    P2T_tic=tic; %phase 2 timer start    
    
    %STEEPEST/GRADIENT DESCENT CODE GOES HERE=================================%
    
    %{
    %Editing + Data Extraction codes at your disposal (use one or the other in your descent code):
    %Subroutine inputs:
   
    schem - already defined at top of pipeline subroutine (inputs{2})
    specs - target frequency defined at top of pipeline subroutine (metrics{1})
    in_param - an array of variable parameter values found in database (candidates(1,2:4))
    param_bounds - doesn't appear to be used at all in the subroutine codes...unnecessary input(?)
    d_param - differential change for the variable parameter Jacobian (dvar)
    %}
    
    %preparing inputs to steepest descent subroutine
    specs        = FT;
    in_param     = candidates(1,2:4);
    param_bounds = [75,75,1.5]; %(not used at all in descent code...)
    d_param      = dvar;
        
    %initiate steepest descent algorithm
    [in_seed,UD]=nesterov_gradient_descent(inputs,specs,in_param,param_bounds,d_param,tol,nit);
    %=========================================================================%
 
    if UD==1 %user is happy with local minimum found by Nesterov...
    
        %if local min has been found, stop pipeflow and indicate solution has been found
        disp(['Evaluating constraints (noise margin, etc)...']);    

        %need to plug in final "in_seed" parameter values to extract data for evaluation

        %update input cell for editing/extraction
        vals{1} = inputs{5}{1};
        vals{2} = inputs{5}{2};
        vals{3} = inputs{5}{3};
        vals{4} = inputs{5}{4};
        vals{5} = inputs{5}{5};
        vals{6} = inputs{5}{6}; 
        vals{7} = in_seed(1); %Wp
        vals{8} = in_seed(2); %Wn
        vals{9} = in_seed(3); %Wn,b
        inputs{5} = vals;%vals

        [raw_data,~]=edit_extract(inputs);

        %prepare inputs to evaluation subroutine for verification
        final_vars=in_seed; %(nm,nm,um)
        final_f0=raw_data.f0; %(MHz)
        final_err=100*abs(final_f0-FT)/FT; %final error (%)

        data_in{1} = raw_data.vout; %evaluate local minimum found
        reqs{1} = SLB;
        reqs{2} = DCB;
        %decide whether discovered point is valid
        [td]=evaluator(data_in,reqs);

        if td==1
            disp(['Local minimum appears valid; send to database...']);
            appendDB=1; %switch on database append
            pipe_switch=4; %append found good design to database
            solution=1; %output good design
        elseif td==0
            disp(['Local minimum is not valid; begin convex module...']);
            localfail=1; %local minimum found is not valid
            pipe_switch=3; %begin data fitting/convex optimization
        end

        P2T=toc(P2T_tic);
        disp([' ']);
        disp(['Phase 2 complete. Time = ',num2str(P2T),' seconds.']);
    
    
    elseif UD==0 %user is not happy with local minimum found by Nesterov...
    
        pipe_switch=0; %stop flow
        solution=0; %no solution found
    
    end
    
end


if pipe_switch==3 %convex optimization

%all the MAIN subroutines go here!!
disp(['----------------------------------------------------------------']);
if localfail==0
    disp(['PHASE TWO: DATA FITTING & CONVEX OPTIMIZATION']);
elseif localfail==1
    disp(['PHASE THREE: DATA FITTING & CONVEX OPTIMIZATION']);
end
P3T_tic=tic; %phase 3 timer start

%DATA FITTING/CONVEXIFICATION/GEOMETRIC PROGRAMMING CODE GOES HERE========%

disp(['WARNING! Convexification module not completed...']);
solution=0; %no solution found

%}
%=========================================================================%

P3T=toc(P3T_tic);
disp(' ');
if localfail==0
    disp(['Phase 2 complete. Time = ',num2str(P3T),' seconds.']);
elseif localfail==1
    disp(['Phase 3 complete. Time = ',num2str(P3T),' seconds.']);
end

end

if pipe_switch==4 %display results and append database
    
    disp(['----------------------------------------------------------------']);
    disp(['FINAL PHASE: APPEND SOLUTION TO DATABASE AND/OR OUTPUT DESIGN']);
    P4T_tic=tic; %phase 4 timer start

%     relerr=100*abs(candidates(1,1)-FT)/FT;
    
    disp([' ']);
    disp(['Good design solution metric:']);
    disp(['f0 = ',num2str(final_f0),' MHz']);
    disp(['Relative error = ',num2str(final_err),'%']);
    disp(['Good design solution parameters:']);
    %will want to automate this using the header in the future!!
    disp(['Wp = ',num2str(final_vars(1)),' nm']);
    disp(['Wn = ',num2str(final_vars(2)),' nm']);
    disp(['Wn,b = ',num2str(final_vars(3)),' um']);  

    disp([' ']);
    if appendDB==1
        disp(['Updating database...']);
        database=importdata(['database.mat']);
        database=[database;[final_f0,final_vars]];
        save(['database.mat'],'database');
    end
    fclose('all'); %close all running instances of LTspice

    P4T=toc(P4T_tic);
    disp([' ']);
    disp(['Phase 4 complete. Time = ',num2str(P4T),' seconds.']);
    
end

fclose('all'); %close all running instances of LTspice

disp(['----------------------------------------------------------------']);
TPT=toc(TPT_tic);
disp(['The total pipeline time (TPT) = ',num2str(TPT),' seconds.']);
disp(' ');

%for troubleshooting purposes (comment out when not in use)
% solution=1;

if solution==1
    disp(['Output found good design...']);
    
    %use winning candidate's parameters (again, will automate this in the
    %future...)
    
    %update input cell for editing/extraction
    vals{1} = inputs{5}{1};
    vals{2} = inputs{5}{2};
    vals{3} = inputs{5}{3};
    vals{4} = inputs{5}{4};
    vals{5} = inputs{5}{5};
    vals{6} = inputs{5}{6}; 
    vals{7} = final_vars(1); %Wp
    vals{8} = final_vars(2); %Wn
    vals{9} = final_vars(3); %Wn,b
    
    %editor inputs
    final_inputs{1} = folder; %folder
    final_inputs{2} = schem; %schem
    final_inputs{3} = comps;%comps
    final_inputs{4} = prefixes;%prefixes
    final_inputs{5} = vals;%vals
    final_inputs{6} = iter;%iter
    
    %generate good design schematic (use good design parameters)
    scpd_editor(final_inputs);    
    
    exe_path=strtrim(exe_path);    
    
    schem_edit=[schem,'_edit'];
    schemfile=[folder,'/',schem_edit,'.asc']; %schematic filename
    % netfile=[folder,'/',schem_edit,'.net']; %corresponding netlist filename
    rawfile=[folder,'/',schem_edit,'.raw']; %raw filename generated by LTspice
    
    %execution commands
    % run_schem=['"' exe_path '" -ascii -b -run "' netfile '"'];
    run_schem=['"' exe_path '" -run "' schemfile '"'];
    system(run_schem);

elseif solution==0
    
    disp(['No solution found, ending program :(']);

end


%option to save all data generated in pipeline function
if savedata==1
    save('pipeline_data.mat');
end


end