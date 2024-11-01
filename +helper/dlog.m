function dlog(msg, source)
    if nargin < 2 % If source is empty/not provided
        source = dbstack();
        source = source(2).name; % Print parent function
    end
    % Put all internal strings through this function, for troubleshooting
    fprintf('[%s]\t%s\n', source, msg);
end

