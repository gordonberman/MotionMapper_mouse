function [pixels,thetas,means,stDevs,vidObjs] = findRadonPixels(filePath,numToTest,parameters)
%findRadonPixels finds a set of pixels that contain much of the observed
%variance in the images given a set of aligned .tiff files.  A gui will
%pop-up asking the user to select a cut-off point.  The left-most minimum
%of the distribution is usually optimal for this
%
%   Input variables:
%
%       filePath -> directory containing aligned .tiff files
%       numToTest -> number of images from which to calculate values
%       parameters -> struct containing non-default choices for parameters
%
%
%   Output variables:
%
%       pixels -> list of the chosen high-variance pixels
%       thetas -> list of the thetas used in the Radon transform
%       means -> double array containing the mean radon-transform values
%       stDevs -> double array containing the standard deviation of the
%                   radon-transform values
%       vidObjs -> VideoReader objects for each of the aligned avi files
%
% (C) Gordon J. Berman, 2016
%     Emory University

    numPixels = 7000;

    addpath(genpath('./utilities/'));
    addpath(genpath('./PCA/'));
    
    if nargin < 3
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
    
    numThetas = parameters.num_Radon_Thetas;
    spacing = 180/numThetas;
    thetas = linspace(0,180-spacing,numThetas);
    scale = parameters.rescaleSize;
    
    
    
    [means,stDevs,vidObjs] = findImageSubsetStatistics(filePath,numToTest,thetas,scale);
    
    
    [Y,X] = hist(stDevs(:),100);
    figure
    %test = true;
    
    semilogy(X,Y,'o-')
    %title('Select Truncation Point (Best if at Local Minimum)','fontsize',16,'fontweight','bold');
    xlabel('Standard Deviation','fontsize',14,'fontweight','bold')
    ylabel('Number of Counts','fontsize',14,'fontweight','bold')
    %     while test
    %
    %         [x,~] = ginput(1);
    %         pixels = find(stDevs(:) > x);
    %         title(['# of Pixels = ' num2str(length(pixels)) '.  Is this OK? [y/n]'],'fontsize',16,'fontweight','bold');
    %
    %         [~,~,button] = ginput(1);
    %         while isempty(button)
    %             [~,~,button] = ginput(1);
    %         end
    %
    %         while button ~= 121 && button ~= 110 && button ~= 89 && button ~= 78
    %             [~,~,button] = ginput(1);
    %         end
    %
    %         if button == 121 || button == 89
    %             test = false;
    %         end
    %
    %     end
    
    vals = sort(stDevs(:),'descend');
    pixels = find(stDevs > vals(numPixels));
    
    if ~isempty(p) && parameters.closeMatPool
        delete(p);
    end

    