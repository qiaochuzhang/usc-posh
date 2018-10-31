function [iter]=scpd_editor(inputs)
%Sub-Circuit Parameter-Directive editor
%edit all relevant schematics

%{
%function testing grounds
clear all; clc; clf; close all;
load('testdata.mat');
%}

% %{
folder   = inputs{1};
schem    = inputs{2};
comps    = inputs{3};
prefixes = inputs{4};
vals     = inputs{5};
iter     = inputs{6};
%}

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

end