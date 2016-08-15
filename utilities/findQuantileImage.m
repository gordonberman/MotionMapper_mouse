function [medianImage,images] = findQuantileImage(vidObj,N,quantileValue,firstFrame,lastFrame)

    if nargin < 3 || isempty(quantileValue)
        quantileValue = .25;
    end

    N = min(N,vidObj.NumberOfFrames);
    if nargin < 4 || isempty(firstFrame)
        firstFrame = 1;
    end
    
    if nargin < 5 || isempty(lastFrame)
        lastFrame = vidObj.NumberOfFrames;
    end

    if vidObj.NumberOfFrames <= N
        
        images = read(vidObj,[1 vidObj.NumberOfFrames]);
    
    else
        
        q = randi(lastFrame - firstFrame + 1) + firstFrame - 1;
        
        image = read(vidObj,q);
        s = size(image);
    
        images = zeros(s(1),s(2),N);
        images(:,:,1) = image(:,:,1);
        
        for i=2:N
            image = read(vidObj,randi(vidObj.NumberOfFrames));
            images(:,:,i) = image(:,:,1);
        end
        
    end
    
    
    
    medianImage = quantile(images,quantileValue,3);