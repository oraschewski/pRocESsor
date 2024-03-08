function [yIndFl,yIndCe,ratioFl,ratioCe,yFlInRange,yCeInRange] = getFloorCeil(y, yMax)
    %GETFLOORCEIL find indices above and belove for an index array that
    %does not only contain integers.
    %   Detailed explanation goes here

    % Find indices above and below y
    yIndFl = floor(y);
    yIndCe = ceil(y);
    
    % Compute distance to y, needed for scaling
    ratioFl = y - yIndFl;
    ratioCe = yIndCe - y;
    
    % Ensure to consider point once if its values is an integer, yCe = yFl and
    % ratioFl = ratioCe = 0
    ratioSum = ratioFl + ratioCe;
    ratioCe(ratioSum == 0) = 1; 
    
    %Find values in range
    yFlInRange = yIndFl > 0 & yIndFl <= yMax;
    yCeInRange = yIndCe > 0 & yIndCe <= yMax;
    
    % Cut values and rations to range
    yIndFl  = yIndFl(yFlInRange);
    ratioFl = ratioFl(yFlInRange);
    yIndCe  = yIndCe(yCeInRange);
    ratioCe = ratioCe(yCeInRange);
end

