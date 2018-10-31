function [grad_t] = grad(inputs,in_param, f_in_param, f_s,d_param)
w=in_param;
f_w = f_in_param;
for i=1:size(d_param,2)
    w_temp = w;
    w_temp(i) = w_temp(i)+d_param(i);
    %call edit extract to get get f(w_temp)
    
    %update input cell for editing/extraction
    vals{1} = inputs{5}{1};vals{2} = inputs{5}{2};vals{3} = inputs{5}{3};
    vals{4} = inputs{5}{4};vals{5} = inputs{5}{5};vals{6} = inputs{5}{6}; 
    vals{7} = w_temp(1); %Wp
    vals{8} = w_temp(2); %Wn
    vals{9} = w_temp(3); %Wn,b
    inputs{5}=vals;
    [raw_data,~]=edit_extract(inputs);    
    
%     f(i) = f(w_temp);
    f(i) = raw_data.f0;
    diff(i) = f(i)-f_w;
end

grad_t = 2*(f_w - f_s)*(diff./d_param);