function [trainingSetData,trainingSetAmps,projectionFiles] = ...
            runEmbeddingSubSampling(projectionDirectory,parameters)
%runEmbeddingSubSampling generates a training set given a set of .mat files
%
%   Input variables:
%
%       projectionDirectory -> directory path containing .mat projection 
%                               files.  Each of these files should contain
%                               an N x pcaModes variable, 'projections'
%       parameters -> struct containing non-default choices for parameters
%
%
%   Output variables:
%
%       trainingSetData -> normalized wavelet training set 
%                           (N x (pcaModes*numPeriods) )
%       trainingSetAmps -> Nx1 array of training set wavelet amplitudes
%       projectionFiles -> cell array of files in 'projectionDirectory'
%
%
% (C) Gordon J. Berman, 2016
%     Emory University
    
    addpath(genpath('./utilities/'));
    addpath(genpath('./t_sne/'));
    
    if nargin < 2
        parameters = [];
    end
    parameters = setRunParameters(parameters);
    
    
    numProcessors = parameters.numProcessors;
    p = gcp('nocreate');
    c = parcluster;
    numAvailableProcessors = c.NumWorkers;
    
    if numProcessors > 1 && isempty(p)
     
        if numAvailableProcessors > numProcessors
            numProcessors = numAvailableProcessors;
            parameters.numProcessors = numAvailableProcessors;
        end
        
        if numProcessors > 1
            p = parpool(numProcessors);
        end
        
        
    else
        
        if numProcessors > 1
            currentNumProcessors = p.NumWorkers;
            numProcessors = min([numProcessors,numAvailableProcessors]);
            if numProcessors ~= currentNumProcessors
                delete(p);
                p = parpool(numProcessors); 
            end
        end
        
    end
    
    
    projectionFiles = findAllImagesInFolders(projectionDirectory,'.mat');
    
    N = parameters.trainingSetSize;
    L = length(projectionFiles);
    numPerDataSet = round(N/L);
    numModes = parameters.pcaModes;
    numPeriods = parameters.numPeriods;
     
    trainingSetData = zeros(numPerDataSet*L,numModes*numPeriods);
    trainingSetAmps = zeros(numPerDataSet*L,1);
    
    for i=1:L
        
        fprintf(1,['Finding training set contributions from data set #' ...
            num2str(i) '\n']);
        
        currentIdx = (1:numPerDataSet) + (i-1)*numPerDataSet;
        
        [yData,signalData,signalAmps,~] = ...
                file_embeddingSubSampling(projectionFiles{i},parameters);
            
        [trainingSetData(currentIdx,:),trainingSetAmps(currentIdx)] = ...
            findTemplatesFromData(signalData,yData,signalAmps,...
                                numPerDataSet,parameters);
            
            
    end
    
    
    if ~isempty(p) && parameters.closeMatPool
        delete(p);
    end
