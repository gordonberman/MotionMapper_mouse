function medianImage = findMedianImage(vidObj,N,firstFrame,lastFrame)

    N = min(N,vidObj.NumberOfFrames);
    if nargin < 3 || isempty(firstFrame)
        firstFrame = 1;
    end
    
    if nargin < 4 || isempty(lastFrame)
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
    
    medianImage = median(images,3);