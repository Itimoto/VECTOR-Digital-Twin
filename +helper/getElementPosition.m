% Snippet from `helperViewArray`
function [elPos, elSpacing] = getElementPosition(cfgArray)
    elSpacing = repmat(cfgArray.ElementSpacing, ...
        [1, 2/length(cfgArray.ElementSpacing)]);
    
    % Calculate element positions along rows and columns
    numCol = cfgArray.Size(2);
    numRow = cfgArray.Size(1);    
    rowPos = (0:numCol-1)*elSpacing(2);
    rowPos = rowPos - rowPos(end)/2;
    colPos = (0:numRow-1)*elSpacing(1);
    if numCol > 1
        colPos =  colPos(end)/2 - colPos;
    else
        colPos =  colPos - colPos(end)/2;
    end
    
    % Formulate the position grid on the plane where the array panel lies
    expRowPos = kron(rowPos, ones(1, numRow));
    expColPos = repmat(colPos, 1, numCol);
    
    % Formulate [x;y;z] positions
    numEl = prod(cfgArray.Size);
    
    elPos = [zeros(1, numEl); expRowPos; expColPos];
end

