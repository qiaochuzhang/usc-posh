function [database]=database_gen(metric_grid,var_arrays)
%database generator

%{
clear all;clc;clf;close all;
var_arrays=importdata(['var_arrays.mat']);
metric_grid=importdata(['metric_grid.mat']);
%}

nvars=size(var_arrays,2);
for i=1:nvars
    len(i)=numel(var_arrays{i});
end

dblen=numel(metric_grid); %total number of valid solutions

%generate initial set of coordinates ("known good solutions")

%preallocation should only be done if number of NaN's are known...
% database=zeros(dblen,1+nvars); %(1 metric + 3 variables);

r=1;
for k=1:len(3)
    for j=1:len(2)
        for i=1:len(1)
            if double(isnan(metric_grid(i,j,k)))==0 %must NOT include failed values
                database(r,:)=[metric_grid(i,j,k),var_arrays{1}(i),var_arrays{2}(j),var_arrays{3}(k)];
                r=r+1;
            end
        end
    end
end

%export initialized database
save(['database.mat'],'database');

end