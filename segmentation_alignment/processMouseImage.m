function [outImage,centroid,upper_left_corner,threshold] = processMouseImage(...
                image,backgroundImage,imageLength,threshold,imageThreshold,...
                dilateSize,fmin,fmax,openSize,isAlbino)
%processMouseImage takes and image, finds a region containing a mouse, and
%           segments the mouse from the background
%
%
% (C) Gordon J. Berman, 2016
%     Emory University
           
    addpath(genpath('../utilities'));        
            
    minThreshold = 5;
    mirrorRange = 40;
    s = size(image);
    backgroundImage = uint8(backgroundImage);
    image = uint8(image);
    
    if nargin < 10 || isempty(isAlbino)
        isAlbino = false;
    end
    
    
    if nargin < 3 || isempty(imageLength)
        imageLength = 201;
    end
    
    
    if nargin < 4 || isempty(threshold)
        if ~isAlbino
            threshold = 50;
        else
            threshold = 5;
        end
    end
    
    
    if nargin < 5 || isempty(imageThreshold)
        if ~isAlbino
            imageThreshold = 150;
        else
            imageThreshold = 1;
        end
    end
    
    if nargin < 6 || isempty(dilateSize)
        dilateSize = 3;
    end
    
    if nargin < 7 
        fmin = [];
    end
    
    if nargin < 8 
        fmax = [];
    end
    
    if nargin < 9 
        openSize = [];
    end
    
    if isAlbino
        z = round(s(1)/2) + (-mirrorRange:mirrorRange);
        w = mean(imcomplement(backgroundImage(z,:)));
        wallThreshold = median(w(w > threshold));
    end
    

    q = 0;
    while sum(q(:)) == 0 && threshold >= minThreshold
        mask = double(backgroundImage) - double(image) > threshold;
        q = imcomplement(uint8(image));
        if isAlbino
            mask = mask & imcomplement(image) < wallThreshold;
            mask = imfill(mask,'holes');
        end    
        q(~mask) = 0;
        
        mask2 = q > imageThreshold;
        mask2 = imdilate(mask2,strel('square',dilateSize));
        mask2 = imfill(mask2,'holes');
        mask2(~mask) = 0;
        q = immultiply(imcomplement(uint8(image)),mask2);
            
        if sum(q(:)) == 0
            threshold = threshold - 5;
        end
        
    end
    
    
    
    
    
    L = floor(imageLength/2);
    imageLength = 2*L + 1;
    

    if sum(q(:)) > 0
        
        props = regionprops(mask2,q,'Area','WeightedCentroid','PixelIdxList');
        
        idx = argmax([props.Area]);
        if length(idx) > 1
            idx = idx(randi(lenght(idx)));
        end
        
        centroid = props(idx).WeightedCentroid;
        c = round(centroid);
        
        newImage = uint8(zeros(size(q)));
        newImage(props(idx).PixelIdxList) = q(props(idx).PixelIdxList);
        
%         if isAlbino
%             a = newImage;
%             mask = a > 0 & a <= wallThreshold;
%             mask = imdilate(mask,strel('disk',dilateSize));
%             CC = bwconncomp(mask);
%             if CC.NumObjects > 0
%                 lengths = returnCellLengths(CC.PixelIdxList);
%                 b = argmax(lengths);
%                 mask2 = false(size(mask));
%                 mask2(CC.PixelIdxList{b}) = true;
%                 mask2 = imfill(mask2,'holes');
%                 newImage = immultiply(mask2,a);
%                 props = regionprops(mask2,newImage,'WeightedCentroid');
%                 centroid = props(1).WeightedCentroid;
%                 c = round(centroid);
%             end
%         end
        
        
        
        iRange = c(2) + (-L:L);
        jRange = c(1) + (-L:L);
        
        if iRange(1) < 1
            iRange = iRange - iRange(1) + 1;
        end
        
        if iRange(end) > s(1)
            iRange = iRange - (iRange(end) - s(1));
        end
        
        if jRange(1) < 1
            jRange = jRange - jRange(1) + 1;
        end
        
        if jRange(end) > s(2)
            jRange = jRange - (jRange(end) - s(2));
        end
        
        upper_left_corner = [iRange(1) jRange(1)];
        
        outImage = newImage(iRange,jRange);
        
        if ~isempty(fmin) && ~isempty(fmax) && ~isempty(openSize) 
            tempMask = outImage > fmin & outImage < fmax;
            tempMask = imopen(tempMask,strel('square',openSize));
            mask = outImage > 0;
            mask(tempMask) = false;
            mask = imdilate(mask,strel('disk',dilateSize));
            mask = returnLargestConnectedComponentImage(mask);
            mask = imfill(mask,'holes');
            outImage(~mask) = 0;
        end
        
    else
        
        outImage = uint8(zeros(imageLength));
        centroid = [-1 -1];
        upper_left_corner = [-1 -1];
        
    end