function cropLunarGeoTIFF(inputFilename, outputFilename, cropFactor)
    % cropLunarGeoTIFF - Crop a lunar GeoTIFF image by a specified factor.
    %
    % Syntax: cropLunarGeoTIFF(inputFilename, outputFilename, cropFactor)
    %
    % Inputs:
    %   inputFilename  - Path to the input lunar GeoTIFF file.
    %   outputFilename - Path to save the cropped GeoTIFF file.
    %   cropFactor     - Factor by which to crop the image (e.g., 0.5 for 50%).
    
    % Step 1: Read the lunar GeoTIFF
    [Z, R] = readgeoraster(inputFilename); % Read the input GeoTIFF

    % Display the original image
    figure;
    imshow(Z, []);
    title('Original Lunar GeoTIFF');

    % Step 2: Get the dimensions of the original image
    [rows, cols, ~] = size(Z);
    
    % Calculate new dimensions based on cropFactor
    newRows = round(rows * cropFactor);
    newCols = round(cols * cropFactor);

    % Step 3: Define cropping limits in pixel coordinates
    rowCenter = round(rows / 2);
    colCenter = round(cols / 2);
    
    rowStart = max(1, rowCenter - floor(newRows / 2));
    rowEnd = min(rows, rowCenter + floor(newRows / 2) - 1);
    colStart = max(1, colCenter - floor(newCols / 2));
    colEnd = min(cols, colCenter + floor(newCols / 2) - 1);

    % Step 4: Crop the image using mapcrop
    croppedImage = mapcrop(Z, R, [rowStart rowEnd], [colStart colEnd]);

    % Step 5: Create new spatial referencing object for the cropped image
    newR = maprefcells(size(croppedImage), ...
                       R.CellExtentInWorldX, ...
                       R.CellExtentInWorldY, ...
                       R.XWorldLimits, ...
                       R.YWorldLimits);
    
    % Update world limits for the cropped area
    newR.XWorldLimits = [R.XWorldLimits(1) + (colStart - 1) * R.CellExtentInWorldX, ...
                         R.XWorldLimits(1) + colEnd * R.CellExtentInWorldX];
    newR.YWorldLimits = [R.YWorldLimits(1) + (rowStart - 1) * R.CellExtentInWorldY, ...
                         R.YWorldLimits(1) + rowEnd * R.CellExtentInWorldY];

    % Step 6: Write the cropped image to a new GeoTIFF
    geotiffwrite(outputFilename, croppedImage, newR);

    % Display the cropped image
    figure;
    imshow(croppedImage, []);
    title('Cropped Lunar GeoTIFF');
end
