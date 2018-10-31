function [in_seed,UD]=nesterov_gradient_descent(inputs,specs,in_param,param_bounds,d_param,tol,nit)
% parambounds is the parameter array bounds for search space
% schem is the input schem cell array
%Set initial conditions
Lmin=0; %no local minimum found before search (by definition)
t = 1;
% tol = 9; %error tolerance (%)
a_0 = 1; 
w_0=in_param; %initial starting parameter

z=(2*rand(1,3)-1).*(param_bounds) + in_param; %random point in neighborhood for Lipschitz constant
%call edit_extract_2 to get frequenct at w_0

%update input cell for editing/extraction
vals{1} = inputs{5}{1};vals{2} = inputs{5}{2};vals{3} = inputs{5}{3};
vals{4} = inputs{5}{4};vals{5} = inputs{5}{5};vals{6} = inputs{5}{6}; 
vals{7} = in_param(1); %Wp
vals{8} = in_param(2); %Wn
vals{9} = in_param(3); %Wn,b
inputs{5}=vals;
[raw_data,~]=edit_extract(inputs);
f_0=raw_data.f0;

%update input cell for editing/extraction
vals{1} = inputs{5}{1};vals{2} = inputs{5}{2};vals{3} = inputs{5}{3};
vals{4} = inputs{5}{4};vals{5} = inputs{5}{5};vals{6} = inputs{5}{6}; 
vals{7} = z(1); %Wp
vals{8} = z(2); %Wn
vals{9} = z(3); %Wn,b
inputs{5}=vals;
[raw_data,~]=edit_extract(inputs);
f_z=raw_data.f0;

f(1)=f_0;

% f_0=f_0(w_0); % frequency at w_0;
%call edit_extract_2 to get frequency at z
% f_z=f_0(z); %frequency at z
f_s = specs; %user specs
l_0= (f_0 - f_s).^2; %this is the cost function to be minimized 
l_z = (f_z - f_s).^2;
%need to define and call grad function 
disp(['Calculating 1 of 2 pre-descent gradients for all variables...']);
grad_w_0 = grad(inputs,w_0,f_0,f_s,d_param);% gradient of cost function at w_0
disp(['Calculating 2 of 2 pre-descent gradients for all variables...']);
grad_z = grad(inputs,z,f_z,f_s,d_param);
epsilon_initial = norm(w_0 - z)./norm(grad_w_0 - grad_z); %initial step size from Lipschitz constraint

%================= Nesterov Accelerated Gradient Descent =================%

t_switch=0;
epsilon(t)=epsilon_initial;
theta{t}=w_0;
w{t}=w_0;
l(t)=l_0;
a(t)=a_0;

disp(['Nesterov iteration, t = ',num2str(t),', f0(',num2str(t),') = ',...
    num2str(f(t)),' MHz (cost L(',num2str(t),') = ',num2str(l(t)),')']);

while(t_switch==0)
    i=0;
    i_switch = 0;
    grad_t{t}=grad(inputs,w{t},f(t),f_s,d_param);
    while(i_switch==0)
        w_temp = w{t} - 2^(-i)*epsilon(t)*grad_t{t};
        %call edit extract on w_temp and get f_temp
        
        %update input cell for editing/extraction
        vals{1} = inputs{5}{1};vals{2} = inputs{5}{2};vals{3} = inputs{5}{3};
        vals{4} = inputs{5}{4};vals{5} = inputs{5}{5};vals{6} = inputs{5}{6}; 
        vals{7} = w_temp(1); %Wp
        vals{8} = w_temp(2); %Wn
        vals{9} = w_temp(3); %Wn,b
        inputs{5}=vals;
        [raw_data,~]=edit_extract(inputs);
        f_temp=raw_data.f0;
        
        l_temp = (f_temp - f_s).^2;
        cond = l(t) - l_temp - 2^(-i)*epsilon(t)*0.5*(norm(grad_t{t}).^2);
        if (cond >= -0.05) %the only heuristic at the moment...
            i_switch=1;
        else
            i=i+1;
            disp(['Condition loop iteration, i = ',num2str(i),', condition = ',num2str(cond)]);
        end
    end

    epsilon(t+1) = 2.^(-i)*epsilon(t);
    theta{t+1}= w{t} - epsilon(t+1)*grad_t{t};
    a(t+1)= 0.5*(1+sqrt(4*a(t)^2+2));
    w{t+1}= theta{t+1} + (((a(t)-1)*(theta{t+1}-theta{t}))./(a(t+1)));
    t=t+1;
    
    %call edit_extract to get f_0(t)    
    %update input cell for editing/extraction
    vals{1} = inputs{5}{1};vals{2} = inputs{5}{2};vals{3} = inputs{5}{3};
    vals{4} = inputs{5}{4};vals{5} = inputs{5}{5};vals{6} = inputs{5}{6}; 
    vals{7} = w{t}(1); %Wp
    vals{8} = w{t}(2); %Wn
    vals{9} = w{t}(3); %Wn,b
    inputs{5}=vals;
    [raw_data,~]=edit_extract(inputs);
    f(t)=raw_data.f0;
    
    %compute cost function (using relative error for a single metric, for now...)
    l(t)=(f(t)-f_s).^2;
%     l(t)=100*abs(f(t)-f_s)/f_s; %units of (%) away from target

    disp(['Nesterov iteration, t = ',num2str(t),', f0(',num2str(t),') = ',...
    num2str(f(t)),' MHz (cost L(',num2str(t),') = ',num2str(l(t)),')']);

    %compute local minimum measure (to test convergence to incorrect minimum)
    if t>nit
        errcon=mean(abs(diff(l((t-nit):t))));
        if errcon<1 && l(t)>tol
            Lmin=1;
            t_switch=1; %break out of while loop
        end
    end
    
    if l(t)<=tol
        t_switch=1; %break out of while loop
    end

end

if t_switch==1
    if Lmin==1
        disp(['Objective unreachable within error tolerance bounds (tolerance = ',...
            num2str(tol),'%; error = ',num2str(l(t)),'%)']);
        UD=input(['Use results from local minimum found? (1=Yes,0=No): ']);
    elseif Lmin==0
        disp(['Convergence to a local minimum detected.']);
        UD=1;
    end
end

%export latest array
in_seed=w{t};


end
 





