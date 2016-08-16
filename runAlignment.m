function outputStruct = runAlignment(fileName,outputPath,startImage,finalImage,parameters)
%runAlignment runs the alignment and segmentation routines on a .avi file
%   and saves the output files to a directorty
%
%   Input variables:
%
%       fileName -> avi file to be analyzed
%       outputPath -> path to which files are saved
%       startImage -> first image in path to be analyzed
%       finalImage -> last image in path to be analyzed
%       parameters -> struct containing non-default choices for parameters
%
%
%   Output variable:
%
%       outputStruct -> struct containing found alignment variables
%
% (C) Gordon J. Berman, 2016
%     Emory University


    addpath(genpath('./utilities/'));
    addpath(genpath('./segmentation_alignment/'));
    
    
    if ~exist(outputPath,'dir')
        mkdir(outputPath);
    end
    
    
    if nargin < 3 || isempty(startImage)
        startImage = 1;
    end
    
    
    if nargin < 4 || isempty(finalImage)
        finalImage = [];
    end
    
    
    if nargin < 5 || isempty(parameters)
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
        
    
    %     if matlabpool('size') ~= parameters.numProcessors;
    %         matlabpool close force
    %         if parameters.numProcessors > 1
    %             matlabpool(parameters.numProcessors);
    %         end
    %     end
    
    
    outputStruct = ...
        alignImages_Radon_parallel_avi(fileName,startImage,finalImage,...
                                        outputPath,parameters);
    
    outputStruct.parameters = parameters;
                                    
                                    
    %See alignImages_Radon_parallel_avi for definitions of these variables                                
    %     outputStruct.Xs = Xs;
    %     outputStruct.Ys = Ys;
    %     outputStruct.angles = angles;
    %     outputStruct.areas = areas;

    %     outputStruct.framesToCheck = framesToCheck;
    %     outputStruct.svdskipped = svdskipped;
    %     outputStruct.areanorm = areanorm;
    
    
    if ~isempty(p) && parameters.closeMatPool
        delete(p);
    end
    
    
    