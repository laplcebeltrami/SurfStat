function [ pval, observed, thresholds ] = SurfStatTransposeT(x,y, per_t);

%Corrected P-values based on the transpostion test.
%
% Usage: [ pval, peak ] = SurfStatTransposeT(x,y,t);
%
% x         = m x v matrix of data, m = #samples in group 1, ; v=#vertices;
% y         = n x v matrix of data, n = #samples in group 1, ; v=#vertices;
%
% per_t       = number of permutations. For sample size below 100 in each
%               group, use per_t =100000. For sample size above 100 in each 
%               group, use per_t = 1000000. 
%
% observed.t      = t-stat value over vertices
% observed.max    = max t-stat value over surface
% observed.min    = min t-stat value over surface
%
% thresholds.max = thresholds corresesponds to the corrected p-value of
% 0.025, 0.05, 0.01 and 0.001 on the right tail
%
% thresholds.min = thresholds corresesponds to the corrected p-value of
% 0.025, 0.05, 0.01 and 0.001 on the left tail
%
% pval.max    = corrected p-value after multiple comparisions correction for right tail
% pval.min    = corrected p-value after multiple comparisions correction for left tail
%
%
% Reference: Chung, M.K., Xie, L, Huang, S.-G., Wang, Y., Yan, J., Shen, L. 
% 2019 Rapid acceleration of the permutation test via transpositions, International 
% Workshop on Connectomics in NeuroImaging, Lecture Notes in Computer Science (LNCS) 
% 11848:42-53. An earlier version with different application

% This code is downloaded from
% https://github.com/laplcebeltrami/SurfStat 
% as the part of SurfStat-Extended
%
% (C) 2022 Moo K. Chung
%  mkchung@wisc.edu
% University of Wisconsin-Madison

%observation
m=size(x,1);
n=size(y,1);

observed.t=(mean(x)-mean(y)).*sqrt(m*n*(m+n-2)./((m+n)*((m-1)*var(x)+(n-1)*var(y))));%observed t-stat
observed.max = max(observed.t);
observed.min = min(observed.t);

%transposition test
per_t=1000;
[stat_t, time_t] = test_transpose_minmax(x,y,per_t); 
%stat_t.max and stat_t.min have max and min t-stat value for all
%transpositions. 

%Treshold correpsoidng to corrected p-values
thresholds.max = quantile(stat_t.max, [0.975 0.95 0.99 0.999]); 
thresholds.min = quantile(stat_t.min, [0.025 0.05 0.01 0.001]);


%----------------------------
% The overall multiple-comparsion corrected p-value (p_max) for the one sided test 
% (right tail) is obtained by solving  observed.max = quantile(stat_t.max, p_max) 

quantitle_max = quantile(stat_t.max, 0:0.00001:1); %5-decimal accuracy
% quantitle plot
% figure; plot( 0:0.00001:1, p_right, 'LineWidth', 2)

ind_max = find(quantitle_max>=observed.max);
if isempty(ind_max)
    pval.max = 0.00001; % that's the accuracy. 
else
    pval.max = quantitle_max(min(ind_max)); %use min of max
end


% The overall multiple-comparsion corrected p-value (p_min) for the one sided 
% test (left tail) is obtained by solving  observed.min = quantile(stat_t.min, p_min) 

quantitle_min = quantile(stat_t.min, 0:0.00001:1); %5-decimal accuracy
% quantitle plot
% figure; plot( 0.5:0.00001:1, p_right, 'LineWidth', 2)

ind_min = find(quantitle_min<=observed.min);

if isempty(ind_max)
    pval.min = 0.00001; % that's the accuracy. 
else
    pval.min = quantitle_min(max(ind_max)); %use max of min
end

%The overall p-value for the two sided test (both tails) is obtained by
% pvalue = p_max + p_min


%----------------------------------------
function [stat, time_t] = test_transpose_minmax(x,y,per_t)
%function [stat_t, time_t] = test_transpose(x,y,per_t)
% The function computes the two-sample t-statitic of vector data using the transposition
% test. If follows the method explained in 
%
% Chung, M.K., Xie, L, Huang, S.-G., Wang, Y., Yan, J., Shen, L. 2019 Rapid 
% acceleration of the permutation test via transpositions, International 
% Workshop on Connectomics in NeuroImaging, in press. 
% http://www.stat.wisc.edu/~mchung/papers/chung.2019.CNI.pdf
%
% INPUT
% x    : input data of size m x l (m= number of subjects, l=number of data per subject) 
% y    : input data of size n x l (n= number of subjects, l=number of data per subject)
% per_t: number of transpositions
%
% OUTPUT
% stat.max:  maximum t-statistic over all features for all transpositions
% stat.min:  minimum t-statistic over all features for all transpositions
% time_t:  run time it took to compute the statistics
%
%
% This code is downloaded from
% http://www.stat.wisc.edu/~mchung/transpositions
%
% (C) 2019 Moo K. Chung, Yixian Wang  
%  mkchung@wisc.edu
% University of Wisconsin-Madison


tic;

m=size(x,1);
n=size(y,1);

l=size(x,2); %dimension of vector
%stat_t=zeros(per_t,l); %computed statistic over transpostions will be saved as stat_t.
stat.max =zeros(per_t,1); 
stat.min =zeros(per_t,1);


pi1=randi(m,1,per_t); %generate random transposition indices in group 1
pi2=randi(n,1,per_t); %generate random transposition indices in group 2

z=[x;y];
z=z(randperm(m+n));
x=z(1:m);y=z(m+1:m+n); %initial random shuffle

%------------------
%initial values of sum and squared sum functions
mx1=sum(x);my1=sum(y);%initial summation (with no division)
mean_x=mean(x); mean_y=mean(y);
vx1=sum((x-mean_x).*(x-mean_x)); vy1=sum((y-mean_y).*(y-mean_y));%initial squared sum


%random sequantical transposition
for i=1:per_t
    
    a=x(pi1(i));x(pi1(i))=y(pi2(i));
    b=y(pi2(i));y(pi2(i))=a;%tranpose one element between x and y
    
    mx2=mx1+b-a;my2=my1+a-b;%update summation function
    
    vx2=vx1+(mx1.^2-mx2.^2)/m+b.^2-a.^2; %update squared sum function
    vy2=vy1+(my1.^2-my2.^2)/n+a.^2-b.^2; 
    
    %stat_t(i,:)=(mx2/m-my2/n).*sqrt(m*n*(m+n-2)./((m+n)*(vx2+vy2)));%update t-statistic  
    stat_t =(mx2/m-my2/n).*sqrt(m*n*(m+n-2)./((m+n)*(vx2+vy2)));%update t-statistic 
    mx1=mx2;my1=my2;vx1=vx2;vy1=vy2;%prepare for the next transposition
    
    stat.max(i) = max(stat_t);
    stat.min(i) = min(stat_t);
    %   i=i+1;
    %   time=toc;
end

time_t = toc;
