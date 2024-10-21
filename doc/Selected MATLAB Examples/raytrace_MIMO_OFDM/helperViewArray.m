function helperViewArray(array)
%HELPERVIEWARRAY Visualize array geometry
%   helperViewArray(ARRAY) visualizes the array geometry of the arrayConfig
%   object, ARRAY.

%   Copyright 2020 The MathWorks, Inc. 

% Get element positions
[elPos, elSpacing] = getElementPosition(array);
arraySize = array.Size;

% Plot elements in y-z plane
figure('Position', [360 360 320 320])
scatter(elPos(2,:), elPos(3,:), 100, 'filled');
hold on; axis equal; grid off;

% Label elements
if arraySize(2) > 1 
    textXOffset = elSpacing(2)/10;
else
    textXOffset = elSpacing(1)/10;
end

if arraySize(1) > 1
    textYOffset = elSpacing(1)/10;
else
    textYOffset = elSpacing(2)/10;
end

for i = 1:size(elPos, 2)
    text(elPos(2,i) + textXOffset, elPos(3,i) - textYOffset, num2str(i));
end

% Remove ticks
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.XTick = [];
ax.XTickLabel = [];
ax.YTick = [];
ax.YTickLabel = [];

title('Array Geometry');
xlabel('y'); ylabel('z');
if arraySize(2) == 1
    xlim([-elSpacing(1), elSpacing(1)]);
else
    xlim([-elSpacing(2)*arraySize(2)/2, elSpacing(2)*arraySize(2)/2]);
end

if arraySize(1) == 1
    ylim([-elSpacing(2), elSpacing(2)]);
else
    ylim([-elSpacing(1)*arraySize(1)/2, elSpacing(1)*arraySize(1)/2]);
end

end

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

% [EOF]