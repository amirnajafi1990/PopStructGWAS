function [P, Z, NumOfSNPs, SNPs, F] = association_mapping(x,y,numOfAlleles,numOfClusters,lambda0,eta,eta_,delta,burnIn,numOfIterations,runNumber)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Main File: Association Mapping algorithm %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% INPUTS: %%%%%%%%%%%%%%%%%%%%%%%
%%% "k" is the number of clusters
%%% "x" is 2d matrix (L*n)
%%% "y" is a 1d matrix (n*1) which is either 1 or 0

%%% "lambda0" is the initial vector of dirichlet parameters for generating P
%%% "eta" is the poisson parameter for generating number of SNPSs
%%% "eta_" represents the number of numOfSNPs to consider
%%% "delta" is the parameter representing the maximum neighborhood distance to look for

%%%%%%%%%%%%%%%%%%% OUTPUTS: %%%%%%%%%%%%%%%%%%%%%%%
%%% "P" is the output allele frequencies
%%% "Z" is the output 

%***************  initializing matrices  *******************
dataSize = size(x);
loci = dataSize(1);
individuals = dataSize(2);

lambda = zeros(loci,numOfAlleles,numOfClusters);
p = zeros(loci,numOfAlleles,numOfClusters);
finalP = zeros(loci,numOfAlleles,numOfClusters);

z=zeros(individuals,1);
z_previous = zeros(individuals,1); %z for previous iteration is checked to speed the code up (p is not changed if z == z_previous)
z_previous_model = zeros(individuals,1); %z for previous time the model was used is checked to speed the code up (model is not changed if z == z_previous_model)
finalZ = zeros(individuals,1);

bestSNPs = -1*ones(numOfClusters,eta+eta_); % Most probable set of SNPs
bestEntropy = -1000000*ones(numOfClusters,1);  % entropy associated to the most probable set of SNPs
bestModelTable = -1*ones(numOfAlleles^(eta+eta_),2,numOfClusters);
bestNumOfSNPs = -1*ones(numOfClusters,1);
bestIndexes = -1*ones(individuals,numOfClusters);
finalNumOfSNPs = zeros(numOfClusters,1);
finalSNPs = zeros(numOfClusters,eta+eta_);
finalF = zeros(numOfAlleles^(eta+eta_),numOfClusters);

modelUsed=1;

%***************  end of initialization  *******************

% ************** generating random z **************
z=ceil(rand(individuals,1)*numOfClusters);
%***************  end of generating z  *******************

% ************** starting iterations **************

for iteration=1:numOfIterations
    if(mod(iteration,50) ==0)
        iteration
    end
    
    
    %*********   sampling p from p|x,y,m,z    ***************
    p = zeros(loci,numOfAlleles,numOfClusters);
    if(~isequal(z,z_previous))
        lambda = repmat(lambda0,[loci 1 numOfClusters]);
        for l = 1:loci
            for i = 1:individuals
                lambda(l,x(l,i),z(i)) = lambda(l,x(l,i),z(i))+1;
            end
        end
    end
    
    for k = 1:numOfClusters
        for l = 1:loci
            p(l,:,k) = dirichletrnd(lambda(l,:,k)); %sampling for each locus
        end
    end
    %*********  end of sampling p from p|x,y,m,z    ***************
    
    
    %*********  sampling m from m|x,y,p,z    ***************
    
    if (isequal(z,z_previous_model)) %z is the same as when the model was used
        modelUsed =1;
        fprintf(' model for previous iteration in run number %d iteration %d\n',runNumber,iteration);
    else
        
        % checking whether every cluster has both cases and controls
        case_controls =zeros(numOfClusters,2);%1 represents controls & 2 represents cases
        for i =1:individuals
            case_controls(z(i),y(i)+1)=case_controls(z(i),y(i)+1)+1;
        end
        
        if(ismember(0,case_controls))
            modelUsed =0;
            fprintf('model not used for run number %d iteration %d\n',runNumber,iteration);
        else
            modelUsed =1;
            
            bestEntropy = -1000000*ones(numOfClusters,1);
            bestSNPs = -1*ones(numOfClusters,eta+eta_);
            bestNumOfSNPs = -1*ones(numOfClusters,1);
            bestModelTable = -1*ones(numOfAlleles^(eta+eta_),2,numOfClusters);
            bestIndexes = -1*ones(individuals,numOfClusters);
            
            numOfSNPs_min=max(1,eta-eta_);
            numOfSNPs_max=eta+eta_;
            
            % searches for SNPs = max(1,eta-eta_): eta+eta_ because these have the
            % highest poisson probability
            
            for numOfSNPs = numOfSNPs_min:numOfSNPs_max
                
                maxPoiss = log10(poisspdf(numOfSNPs,eta));
                
                if(maxPoiss>min(bestEntropy))% for speeding up
                    
                    fixedCombinations = combinator(delta-1,numOfSNPs-1,'c'); %fixed combinations given this number of SNPs and delta
                    temp = size(fixedCombinations);
                    numOfSets = temp(1); % num of sets of SNPs which have numOfSNPs elements
                   
                    breakFlag = 0;
                    
                    for l = 1:loci-max(numOfSNPs,1)+1 % investigating different combinations given numOfSNPs and delta
                        if(breakFlag ==1)
                            break;
                        end
                        for s = 1:max(numOfSets,1)
                            
                            if(maxPoiss <= min(bestEntropy))% for speeding up
                                breakFlag =1; % breaks out of the first loop
                                break; % breaks out of the second loop
                            end
                            
                            % performing computatons for finding the best model
                            if(numOfSets ==0)
                                sampleCombination = l;
                            else
                                sampleCombination = [l fixedCombinations(s,:)+l];
                            end
                            
                            if (sampleCombination(numOfSNPs) <= loci) %assumes that combinator sorts the data in ascending manner
                                
                                modelTable=zeros((numOfAlleles^numOfSNPs),2,numOfClusters);
                                indexes = zeros(individuals,1); % inicates what combination of this set each person has
                                
                                for i=1:individuals
                                    radix = 0;
                                    for c=1:numOfSNPs
                                        radix=(x(sampleCombination(c),i)-1)*(10^(numOfSNPs-c))+radix;
                                    end
                                    indexes(i,:) = base2dec(num2str(radix),numOfAlleles)+1; % index of this combination in modelTable for individual i
                                    modelTable(indexes(i),y(i)+1,z(i))=modelTable(indexes(i),y(i)+1,z(i))+1;
                                end
                                
                                mult = log10(poisspdf(numOfSNPs,eta)*ones(numOfClusters,1));
                                
                                for i =1:individuals
                                    mult(z(i)) = mult(z(i))+log10(modelTable(indexes(i),y(i)+1,z(i))/(modelTable(indexes(i),1,z(i)) + modelTable(indexes(i),2,z(i))));
                                end
                                
                                for k=1:numOfClusters
                                    if(bestEntropy(k)<mult(k))
                                        bestEntropy(k)=mult(k);
                                        bestSNPs(k,:) = -1*ones(1,eta+eta_);
                                        bestSNPs(k,1:numOfSNPs)=sampleCombination;
                                        bestNumOfSNPs(k) = numOfSNPs;
                                        bestModelTable(:,:,k)=-1*ones(numOfAlleles^(eta+eta_),2);
                                        bestModelTable(1:numOfAlleles^numOfSNPs,:,k) = modelTable(:,:,k);
                                        bestIndexes(:,k) = indexes;
                                    end
                                end
                                
                            end
                            
                        end
                    end
                end
            end
        end
    end
    
    %*********  end of sampling m from m|x,y,p,z    ***************
    
    %*********** sampling z from z|x,p *************
    
    logProbs = -1000000*ones(individuals,numOfClusters); % log of prob(z(i)=k|x,p) (i*k)
    cumulativeProbs = zeros(individuals,numOfClusters);
    probs = -1*ones(individuals,numOfClusters);          %prob(z(i)=k|x,p) (i*k)
    
    for i = 1:individuals
        z_previous(i)=z(i);
    end
    
    if(modelUsed==1)
        for i = 1:individuals
            z_previous_model(i)=z(i);
        end
    end
    
    for i = 1:individuals
        
        if(modelUsed ==1)
            modelUsedForThisInd =ones(1,numOfClusters); % model is used for this individual (the probability for assigning to each cluster can be computer)
            for k=1:numOfClusters
                
                if((bestModelTable(bestIndexes(i,k),1,k) + bestModelTable(bestIndexes(i,k),2,k)) == 0)
                    modelUsedForThisInd(1,k) =0;
                end
            end
        else
            modelUsedForThisInd =zeros(1,numOfClusters);
        end
        
        if(modelUsedForThisInd==ones(1,numOfClusters))
            
            for k=1:numOfClusters
                mult = log10(bestModelTable(bestIndexes(i,k),y(i)+1,k)/(bestModelTable(bestIndexes(i,k),1,k) + bestModelTable(bestIndexes(i,k),2,k)));
                for l=1:loci
                    mult= mult+log10(p(l,x(l,i),k));
                end
                logProbs(i,k) = mult;
            end
        elseif(modelUsedForThisInd==zeros(1,numOfClusters))
            for k=1:numOfClusters
                mult = 0;
                for l=1:loci
                    mult= mult+log10(p(l,x(l,i),k));
                end
                logProbs(i,k) = mult;
            end
        else
            for k=1:numOfClusters
                if(modelUsedForThisInd(1,k)==1)
                    mult = log10(bestModelTable(bestIndexes(i,k),y(i)+1,k)/(bestModelTable(bestIndexes(i,k),1,k) + bestModelTable(bestIndexes(i,k),2,k)));
                    for l=1:loci
                        mult= mult+log10(p(l,x(l,i),k));
                    end
                    logProbs(i,k) = mult;
                else
                    logProbs(i,k) = -inf;
                end
            end
        end
        
        for k =1:numOfClusters
            if (logProbs(i,k)== -inf)
                probs(i,k) = 0;
            else
                tmp =0;
                for kk =1:numOfClusters
                    tmp = tmp+10^(logProbs(i,kk)-logProbs(i,k));
                end
                probs(i,k) = 1/tmp;
            end
        end
                
        for k =1:numOfClusters
            cumulativeProbs(i,k) = sum(probs(i,1:k),2)./ sum(probs(i,:),2);
        end
        
        for k =1:numOfClusters
            if(rand <= cumulativeProbs(i,k))
                z(i) = k;
                break;
            end
        end
    end
    
    %*********** end of sampling z from z|x,p *************
    
    
    %*************** finalizing the results *************
    if(iteration > burnIn)
        finalP = (finalP*(iteration-burnIn-1)+p)/(iteration-burnIn);
        finalZ = (finalZ*(iteration-burnIn-1)+z)/(iteration-burnIn);
        
        if(modelUsed ==1)
            finalNumOfSNPs = (finalNumOfSNPs*(iteration-burnIn-1)+bestNumOfSNPs)/(iteration-burnIn);
            finalSNPs = (finalSNPs*(iteration-burnIn-1)+bestSNPs)/(iteration-burnIn);
            
            for k = 1:numOfClusters
                if( ~ismember(0,bestModelTable(:,1,k) + bestModelTable(:,2,k)))
                    finalF(:,k) = ((finalF(:,k)*(iteration-burnIn-1))+(bestModelTable(:,2,k)./(bestModelTable(:,1,k)+bestModelTable(:,2,k))))/(iteration-burnIn);
                else
                    fprintf('F and S not changed because of missing data in run number %d iteration %d\n',runNumber,iteration)
                end
            end
        else
            fprintf('F and S not changed since no model used in run number %d iteration %d\n',runNumber,iteration);
        end
    end
end
%***************  end of markov chain iteration  *******************

end

