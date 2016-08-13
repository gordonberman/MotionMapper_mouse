function medianImage = findMedianImage(vidObj,N,firstFrame)

    N = min(N,vidObj.NumberOfFrames);
    if nargin < 3 || isempty(firstFrame)
        firstFrame = 1;
    end

    if vidObj.NumberOfFrames <= N
        
        images = read(vidObj,[1 vidObj.NumberOfFrames]);
    
    else
        
        q = -1;
        while q < firstFrame
            q = randi(vidObj.NumberOfFrames);
        end
        
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