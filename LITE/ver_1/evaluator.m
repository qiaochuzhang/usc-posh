function [td]=evaluator(data_in,reqs)
%evaluation & database storage function

%data generated simulated
vout = data_in{1};

%requirements
SLB = reqs{1};
% DCB = reqs{2};

%compute required metrics
[SL]=state_levels(vout);
% [DC]=dutycycle(vout); %will properly be incorporated in a future version

%test whether metrics fall within constraint bounds (requirements)
if SL(1)<=SLB(1) && SL(2)>=SLB(2)
    td1=1;
else
    td1=0;
end

% if min(DC)>=DCB(1) && max(DC)<=DCB(2)
    td2=1;
% else
%     td2=0;
% end

td = td1*td2; %test decision

end