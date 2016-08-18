function [Xs,Ys,angles,areas,svdskips,centroids,upper_left_corners] = ...
    align_subroutine_parallel_avi(grouping,initialPhi,segmentationOptions,...
                image_path,readout,processorNum,vidObj,areanorm,initialImage)
%align_subroutine_parallel_avi is a subroutine used within
%alignImages_Radon_parallel_avi in order to parallelize properly
%
%
% (C) Gordon J. Berman, 2016
%     Emory University



warning off MATLAB:audiovideo:aviread:FunctionToBeRemoved;
L = length(grouping);
Xs = zeros(L,1);
Ys = zeros(L,1);
angles = zeros(L,1);
areas = zeros(L,1);
svdskips = zeros(L,1);
centroids = zeros(L,2);
upper_left_corners = zeros(L,2);

%fid = fopen(['/Users/gberman/Desktop/' num2str(processorNum) '.txt'],'w');

dilateSize = segmentationOptions.dilateSize;
imageThreshold = segmentationOptions.imageThreshold;
spacing = segmentationOptions.spacing;
pixelTol = segmentationOptions.pixelTol;
basis = segmentationOptions.referenceImage;
openSize = segmentationOptions.openSize;
backgroundImage = segmentationOptions.backgroundImage;
aboveBackgroundThreshold = segmentationOptions.aboveBackgroundThreshold;
isAlbino = segmentationOptions.isAlbino;
fmin = segmentationOptions.fmin;
fmax = segmentationOptions.fmax;
outputImageSize = segmentationOptions.outputImageSize;
maxNumImagesToLoad = segmentationOptions.maxNumImagesToLoad;

startFrame = 0;
switchFrame = maxNumImagesToLoad + 1;
if grouping(end) - grouping(1) <= maxNumImagesToLoad
    vidChunk = read(vidObj,[grouping(1) grouping(end)]);
else
    q = maxNumImagesToLoad + grouping(1) - 1;
    vidChunk = read(vidObj,[grouping(1) q]);
    %fprintf(fid,'%6i\t %6i\t %6i\n',grouping(1),q,grouping(end));
end

if length(size(vidChunk)) == 4
    vidChunk = squeeze(vidChunk(:,:,1,:));
end

open(image_path);
writeVideo(image_path,initialImage);



s = size(basis);
currentPhi = initialPhi;

for j=2:L
    
    if mod(j,readout) == 0
        fprintf(1,'\t Processor #%2i, Image #%7i of %7i\n',processorNum,j,L);
    end
    
    if j == switchFrame
        
        startIdx = grouping(1) + j - 1;
        endIdx = startIdx + maxNumImagesToLoad - 1;
        endIdx = min([endIdx,grouping(end)]);
        vidChunk = read(vidObj,[startIdx endIdx]);
        %fprintf(fid,'%6i\t %6i\t %6i\n',startIdx,endIdx,grouping(end));
        startFrame = j - 1;
        switchFrame = switchFrame + maxNumImagesToLoad;
        
        if length(size(vidChunk)) == 4
            vidChunk = squeeze(vidChunk(:,:,1,:));
        end
        
    end
    
    originalImage = vidChunk(:,:,j-startFrame);
    sCurrent = size(originalImage);
    
    if sCurrent(1) < s(1) || sCurrent(2) < s(2)
        zz = uint8(zeros(s)+255);
        zz(1:sCurrent(1),1:sCurrent(2)) = originalImage;
        originalImage = zz;
    end
    
    [imageOut,centroids(j,:),upper_left_corners(j,:),~] = ...
        processMouseImage(originalImage,backgroundImage,...
        outputImageSize,aboveBackgroundThreshold,...
        imageThreshold,dilateSize,fmin,fmax,openSize,isAlbino);
    
    
    imageOut = rescaleImage(imageOut,areanorm);
    
    areas(j) = sum(imageOut(:) ~= 0);
    
    if max(imageOut(:)) > 0
        
        
        [angles(j),Xs(j),Ys(j),loopImage] = ...
            alignTwoImages(basis,imageOut,currentPhi,spacing,pixelTol,false,imageOut);
        
        currentPhi = angles(j);
        
        if ~isempty(fmax)
            props = regionprops(loopImage > 0,'Centroid');
            props_small = regionprops(loopImage > 0 & loopImage < fmax,'Centroid');
            
            if ~isempty(props)
                xCentroid = props(1).Centroid;
            else
                xCentroid = [0 0];
                
            end
            if ~isempty(props_small)
                xCentroid_small = props_small(1).Centroid;
            else
                xCentroid_small = [0 0];
            end
            
            if ((xCentroid_small(1) > xCentroid(1)) && xCentroid_small(1) > 0) || xCentroid(1) == 0
                
                testPhi = mod(currentPhi+180,360);
                [tempAngle,tempX,tempY,tempImage] = ...
                    alignTwoImages(basis,imageOut,testPhi,spacing,pixelTol,false,imageOut);
                
                props_new = regionprops(tempImage > 0,'Centroid');
                if ~isempty(props_new)
                    xCentroid_new = props_new(1).Centroid;
                else
                    xCentroid_new = [0 0];
                end
                
                
                if (xCentroid_new(1) < xCentroid(1) && xCentroid_new(1) > 0) ...
                        || (xCentroid_new(1) > 0 && xCentroid(1) == 0)
                    angles(j) = tempAngle;
                    Xs(j) = tempX;
                    Ys(j) = tempY;
                    loopImage = tempImage;
                    currentPhi = angles(j);
                end
                
            end
        end
        
        writeVideo(image_path,loopImage);
        
    else
        
        angles(j) = angles(j-1);
        svdskips(j) = 1;
        loopImage = uint8(zeros(size(imageOut2)));
        writeVideo(image_path,loopImage);
        
    end
    
end

%fclose(fid);
