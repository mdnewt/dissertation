%This file generates folder structures for testing
%Script to be run from 'master' directory
%'master' is prepopulated with pre-processed data
%Corresponding masks and stacks are assumed to have same name
% n = total number of atlases (multiple of 4)

n = 50;
p = 'Analysis\';
mkdir(p);
whichiswhich = cell(n,2); % (:,1) is training, (:,2) is unknown

for x = 1:n    
    %%
    % Create base directories numbered from 1 to n (n = total number of atlases)
    mkdir([p,num2str(x)]);
    ISA_mkdirs([p,num2str(x)]);

%     % ** Randomization (comment out if rando.mat already exists) **
%     % Switch: how many samples to include in atlas? 
%     % assumes 40 atlases made in total and 10 atlases for each of the 4 pool sizes
%     if(x > 40)
%         poolSize = 20;
%     elseif(x > 30)
%         poolSize = 10;
%     elseif(x > 20)
%         poolSize = 5 ;
%     elseif(x > 10)
%         poolSize = 3 ;
%     else
%         poolSize = 1 ;
%     end
% 
%     % Divy samples into training and atlas
%     idx = randperm(50)';
%     % **TO DO - initialize randomizer beforehand so it produces the same 
%     %   sets of random numbers each time
%     whichiswhich(x,:) = {sort(idx(1:poolSize)),sort(idx(poolSize+1:end))};
%     save('rando.mat','whichiswhich');
    
end