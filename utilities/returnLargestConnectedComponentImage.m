function outImage = returnLargestConnectedComponentImage(image)

    CC = bwconncomp(image > 0);
    lengths = returnCellLengths(CC.PixelIdxList);
    idx = argmax(lengths);
    outImage = image;
    outImage(:) = 0;
    outImage(CC.PixelIdxList{idx(1)}) = image(CC.PixelIdxList{idx(1)});