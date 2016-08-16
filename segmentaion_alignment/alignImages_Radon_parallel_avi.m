function outputStruct = ...%[Xs,Ys,angles,areas,parameters,framesToCheck,svdskipped,backgroundImage] = ...
    alignImages_Radon_parallel_avi(file_path,startImage,finalImage,image_path,parameters)
%alignImages_Radon_parallel_avi runs the alignment and segmentation routines on a .avi file
%   and saves the output files to a directorty (called by ../runAlignment.m)
%
%   Input variables:
%
%       file_path -> avi file to be analyzed
%       startImage -> first frame of the avi file to be analyzed
%       finalImage -> last frame of the avi file to be analyzed
%       image_path -> path to which files are saved
%       parameters -> struct containing parameters
%
%
%   Output variables:
%
%       Xs -> alignment x translations
%       Ys -> alignment y translations
%       angles -> alignment rotations
%       areas -> segmented areas after segmentation
%       framesToCheck -> frames where a large rotation occurs in a single 
%                           frame.  This might signal a 180 degree rotation
%                           error
%       svdskipped -> blank frames where alignment is skipped
%       backgroundImage -> background image used in image processing
%
% (C) Gordon J. Berman, 2016
%     Emory University

    warning off MATLAB:polyfit:RepeatedPointsOrRescale;
    warning off MATLAB:audiovideo:aviinfo:FunctionToBeRemoved;
    
    
    readout = 100;
    %nDigits = 8;
    
    spacing = parameters.alignment_angle_spacing;
    pixelTol = parameters.pixelTol;
    %minArea = parameters.minArea;
    %asymThreshold = parameters.asymThreshold;
    %symLine = parameters.symLine;
    %initialPhi = parameters.initialPhi;
    dilateSize = parameters.dilateSize;
    %cannyParameter = parameters.cannyParameter;
    imageThreshold = parameters.imageThreshold;
    %maxAreaDifference = parameters.maxAreaDifference;
    %segmentationOff = false;
    basisImage = imread(parameters.basisImagePath);
    %bodyThreshold = parameters.bodyThreshold;
    numProcessors = parameters.numProcessors;
    %rangeExtension = parameters.rangeExtension;
    backgroundImageQuantile = parameters.backgroundImageQuantile;
    
    isAlbino = parameters.isAlbino;
    medianImageNumber = parameters.medianImageNumber;
    outputImageSize = parameters.outputImageSize;
    aboveBackgroundThreshold = parameters.aboveBackgroundThreshold;
    openSize = parameters.openSize;
    
    
    %imageLength,threshold,imageThreshold,dilateSize,isAlbino
    
    %Choose starting and finishing images
    
    vidObj = VideoReader(file_path);
    nFrames = vidObj.NumberOfFrames;
    
    if isempty(startImage)
        startImage = 1;
    end
    
    
    if isempty(finalImage)
        finalImage = nFrames;
    end
    
    
    outputStruct.file_path = file_path;
    outputStruct.startImage = startImage;
    outputStruct.finalImage = finalImage;
    outputStruct.nFrames = nFrames;
    
    
    fprintf(1,'Finding Median Image\n');
    [backgroundImage,testImages] = findQuantileImage(vidObj,medianImageNumber,...
        backgroundImageQuantile,startImage,finalImage);
    
    
    outputStruct.backgroundImage = backgroundImage;
    
    
    segmentationOptions.imageThreshold = imageThreshold;
    %segmentationOptions.cannyParameter = cannyParameter;
    segmentationOptions.dilateSize = dilateSize;
    %segmentationOptions.minArea = minArea;
    segmentationOptions.spacing = spacing;
    segmentationOptions.pixelTol = pixelTol;
    %segmentationOptions.maxAreaDifference = maxAreaDifference;
    %segmentationOptions.segmentationOff = segmentationOff;
    %segmentationOptions.asymThreshold = asymThreshold;
    %segmentationOptions.symLine = symLine;
    segmentationOptions.outputImageSize = outputImageSize;
    segmentationOptions.backgroundImage = backgroundImage;
    segmentationOptions.isAlbino = isAlbino;
    segmentationOptions.openSize = openSize;
    segmentationOptions.backgroundImage = backgroundImage;
    segmentationOptions.aboveBackgroundThreshold = aboveBackgroundThreshold;
    segmentationOptions.imageLength = imageLength;
    
    
    %Find segmenatation parameters
    fprintf(1,'Finding Segmentation Parameters\n');
    medianImageNumber = length(testImages(1,1,:));
    testSegmentedImages = zeros(outputImageSize,outputImageSize,medianImageNumber);
    for i=1:medianImageNumber
        [testSegmentedImages(:,:,i),~,~,~] = processMouseImage(...
                testImages(:,:,i),backgroundImage,imageLength,...
                aboveBackgroundThreshold,imageThreshold,dilateSize,...
                [],[],[],isAlbino);
    end
    
    
    q = double(testSegmentedImages(testSegmentedImages > 0));
    if length(q) > 50000
        q = q(radnperm(length(q),50000));
    end
    obj = gmixPlot(q,3,[],[],true,false,[],[],3);
    qq = linspace(0,255,10000);
    posts = posterior(obj,qq');
    [~,a] = sort(obj.mu);
    b = a(2);
    posts = posts(:,b);
    f = fit(qq',posts(:,b),'linearinterp');
    f0=qq(find(posts(:,b) > .5 & qq' > obj.mu(a(1)),1,'first'));
    f02=qq(find(posts(:,b) > .5 & qq' < obj.mu(a(3)),1,'last'));
    fmin = fzero(@(x) f(x)-.5,f0);
    fmax = fzero(@(x) f(x)-.5,f02);
    
    segmentationOptions.fmin = fmin;
    segmentationOptions.fmax = fmax;

    outputStruct.fmin = fmin;
    outputStruct.fmax = fmax;
    outputStruct.openSize = openSize;
    
    testAreas = zeros(medianImageNumber,1);
    for i=1:medianImageNumber
        [testSegmentedImages(:,:,i),~,~,~] = processMouseImage(...
            testImages(:,:,i),backgroundImage,imageLength,...
            aboveBackgroundThreshold,imageThreshold,dilateSize,...
            fmin,fmax,openSize,isAlbino);
        testAreas(i) = sum(sum(testSegmentedImages(:,:,i) > fmax));
    end
    
    medianArea = median(testAreas);
    basisArea = sum(basisImage(:) > fmax);
    areanorm = sqrt(1.1*medianArea/basisArea);
    outputStruct.areanorm = areanorm;
    

    
    if isempty(image_path)
        image_path   = input('Image Path = ?:  ', 's');
    end
    
    [status,~]=unix(['ls ' image_path]);
    if status == 1
        unix(['mkdir ' image_path]);
    end
    
    
    referenceImage = basisImage;
    
    
    %     [ii,~] = find(referenceImage > 0);
    %     minRangeValue = min(ii) - rangeExtension;
    %     maxRangeValue = max(ii) + rangeExtension;
    
    segmentationOptions.referenceImage = referenceImage;
    %     segmentationOptions.minRangeValue = minRangeValue;
    %     segmentationOptions.maxRangeValue = maxRangeValue;
    
    
    
    %define groupings
    imageVals = startImage:finalImage;
    numImages = length(imageVals);
    minNumPer = floor(numImages / numProcessors+1e-20);
    remainder = mod(numImages,numProcessors);
    count = 1;
    groupings = cell(numProcessors,1);
    for i=1:numProcessors
        if i <= remainder
            groupings{i} = imageVals(count:(count+minNumPer));
            count = count + minNumPer + 1;
        else
            groupings{i} = imageVals(count:(count+minNumPer-1));
            count = count + minNumPer;
        end
    end

    
    
    % Write Out Grouping Start and Finish indices
    groupidx = zeros(length(groupings),2);
    for i = 1:length(groupings)
        groupidx(i,1) = groupings{i}(1);
        groupidx(i,2) = groupings{i}(end);
    end
    
    
    %initialize new avi files
    alignmentFiles = cell(numProcessors,1);
    fDigits = ceil(log10(numProcessors+1e-10));
    for i=1:numProcessors
        qq = num2str(i);
        qq2 = [repmat('0',1,fDigits - length(qq)) qq];
        alignmentFiles{i} = VideoWriter([image_path '/' qq2 '.avi']);
    end
    
    
    x1s = zeros(numProcessors,1);
    y1s = zeros(numProcessors,1);
    angle1s = zeros(numProcessors,1);
    area1s = zeros(numProcessors,1);
    svdskip1s = zeros(numProcessors,1);   
    currentPhis = zeros(numProcessors,1);
    centroid1s = zeros(numProcessors,2);
    ULC1s = zeros(numProcessors,2);
    
    %initialize First Images
    
    images = cell(numProcessors,1);
    for j=1:numProcessors
        
        i = groupings{j}(1);
        originalImage = read(vidObj,i);
        if length(size(originalImage)) == 3
            originalImage = originalImage(:,:,1);
        end
        
        [imageOut,centroid1s(j,:),ULC1s(j,:)] =  ...
            processMouseImage(originalImage,backgroundImage,...
                        imageLength,aboveBackgroundThreshold,...
                        imageThreshold,dilateSize,fmin,fmax,openSize,isAlbino);                   
                    
        [angle_0,x_0,y_0,aligned_0] = ...
            alignTwoImages(referenceImage,imageOut,0,spacing,pixelTol,false,imageOut);          
                    
                    
        props = regionprops(aligned_0 > 0,'Centroid');
        if ~isempty(props)
            x1 = props(1).Centroid;
        else
            x1 = [-1 -1];
        end

        [angle_180,x_180,y_180,aligned_180] = ...
            alignTwoImages(referenceImage,imageOut,180,spacing,pixelTol,false,imageOut);
        
        props = regionprops(aligned_0 > 0,'Centroid');
        if ~isempty(props)
            x2 = props(1).Centroid;
        else
            x2 = [-1 -1];
        end
        
        if (x1(1) < x2(1) && x1(1) > 0) || (x1(1) > 0 && x2(1) < 0)   
            
            images{j} = aligned_0;
            area1s(j) = sum(aligned_0(:) > 0);
            currentPhis(j) = angle_0;
            svdskip1s(j) = area1s(j) == 0;
            x1s(j) = x_0;
            y1s(j) = y_0;
            angle1s(j) = angle_0;
            
        else
            
            images{j} = aligned_180;
            area1s(j) = sum(aligned_180(:) > 0);
            currentPhis(j) = angle_180;
            svdskip1s(j) = area1s(j) == 0;
            x1s(j) = x_180;
            y1s(j) = y_180;
            angle1s(j) = angle_180;
            
        end
        
        
    end
        
   
    
    
    fprintf(1,'Aligning Images\n');
    
    Xs_temp = cell(numProcessors,1);
    Ys_temp = cell(numProcessors,1);
    Angles_temp = cell(numProcessors,1);
    Areas_temp = cell(numProcessors,1);
    svdskips_temp = cell(numProcessors,1);
    centroids_temp = cell(numProcessors,1);
    upper_left_corners_temp = cell(numProcessors,1);
    
    
    parfor i=1:numProcessors
        
        %[Xs_temp{i},Ys_temp{i},Angles_temp{i},Areas_temp{i},svdskips_temp{i}] = ...
        %    align_subroutine_parallel_avi(groupings{i},currentPhis(i),...
        %    segmentationOptions,nDigits,file_path,alignmentFiles{i},readout,i,...
        %    asymThreshold,area1s(i),vidObj,[],areanorm,images{i});
        
        [Xs_temp{i},Ys_temp{i},Angles_temp{i},Areas_temp{i},...
            svdskips_temp{i},centroids_temp{i},upper_left_corners_temp{i}] = ...
        align_subroutine_parallel_avi(groupings{i},currentPhis(i),segmentationOptions,...
                alignmentFiles{i},readout,i,vidObj,areanorm,images{i});
        
        Xs_temp{i}(1) = x1s(i);
        Ys_temp{i}(1) = y1s(i);
        Areas_temp{i}(1) = area1s(i);
        Angles_temp{i}(1) = angle1s(i);
        svdskips_temp{i}(1) = svdskip1s(i);
        centroids_temp{i}(1,:) = centroid1s(i,:);
        upper_left_corners_temp{i}(1,:) = ULC1s(i,:);
        
        close(alignmentFiles{i});   
        
    end
    
    
    outputStruct.Xs = cell2mat(Xs_temp);
    outputStruct.Ys = cell2mat(Ys_temp);
    outputStruct.angles = cell2mat(Angles_temp);
    outputStruct.areas = cell2mat(Areas_temp);
    outputStruct.svdskips = cell2mat(svdskips_temp);
    outputStruct.upper_left_corners = cell2mat(upper_left_corners_temp);
    outputStruct.centroids = cell2mat(centroids_temp);
    
        
    
    x = abs(diff(unwrap(angles.*pi/180).*180/pi));
    outputStruct.framesToCheck = find(x > 90) + 1;
    outputStruct.svdskipped = find(outputStruct.svdskips == 1);
    

    
    

