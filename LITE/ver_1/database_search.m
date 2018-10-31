function [candidates]=database_search(metrics)
%Sort and search database according to weighted metrics (f0) and distance
%bound (eps). This has the potential to become significantly more
%sophisticated in the future.
    
disp('Searching database for compatible solutions...');

%{
FT = 365; %desired output frequency
eps = 5; %test case
eps = eps/100; %convert to number
%}        

FT  = metrics{1}; %target frequency (MHz)
SLB = metrics{2}; %state level bounds (V)
DCB = metrics{3}; %duty cycle bounds
eps = metrics{4}(1); %percent distance of target frequency to search in database (%)
tol = metrics{4}(2); %percent tolerance for gradient/steepest descent (%)
nit = metrics{4}(3); %number of iterations to scan for convergence to incorrect minimum

%load up the database for searching
DB = importdata(['database.mat']);

%sort database according to prioritized metric
m_ind=1; %metric index (will eventually be determined from header and priorities)
SDB=sortrows(DB,m_ind);

%extract all rows within bounds (sorted in increasing frequency)
presorted = SDB( SDB(:,1)>(1-eps/100)*FT & SDB(:,1)<(1+eps/100)*FT,:); 

%sort in order of closest match to target
temp=abs(presorted(:,1)-FT);
[~,ind]=sort(temp);

%sorted candidates (beginning with closest match)
candidates=presorted(ind,:);

end