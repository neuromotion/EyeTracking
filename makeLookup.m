% Caleb R Tulloss
% Brown Neuromotion Lab
% Eye Orientation Tracker
% Angle Lookup Table Generation Script
% Version 14
% 3/24/17

% 4/15/18 note: at this point in development, alpha and beta
% were used instead of theta and phi
% alpha = -phi and beta = theta

% X in this script corresponds to Z in the thesis
% Y in this script corresponds to X in the thesis
% Z in this script corresponds to Y in the thesis

% Plot labels have been changed to be consistent with thesis convention

function makeLookup
%% Constant Information
% Start angle of resonator on eye
resAlpha = 0.575;
% Vector from origin to receiver
recLoc = [1.8, 0, 0];
% Radius of eyeball
radius = 1;

% Number of entries
p = 2^11;

% Number of alpha and beta values
a = 2^8;
b = 2^8;
% Number of angle combinations
size = a*b;

alphaRange = 0.9;
betaRange = 0.9;

% Step sizes based on hemisphere
alphastep = alphaRange*2/a;
betastep = betaRange*2/b;

% B field components generated for each angle
% Used to find the distribution of possible B-fields
% and thus calculate where divisions should fall
Bx = zeros(size, 1);
By = zeros(size, 1);
Bz = zeros(size, 1);

% Length of list to store angles that produce the same index
listlength = 200;

% Compute divisions for X component
xdivpower = 6;
xnumdivisions = 2^xdivpower;

% Compute divisions for Y component
ynumdivisions = p/xnumdivisions;

% The empty lookup table
% Dim 1: X index
% Dim 2: Y index
% Dim 3: Alpha/Beta
% Dim 4: List to store all angles in bucket
lookup = zeros(xnumdivisions, ynumdivisions, 2, listlength);

%% Distribution Determination

% 4/15/18 note: looking back I really should have done these
% calculations with matrix operations instead of loops...

% X component distribution-finding loop
for i = 1:a
    % Calculate alpha in radians
    alpha = resAlpha - alphaRange +  i*alphastep;
    % X component of coil orientation
    Wx = sin(alpha);
    
    for j = 1:b
        % Calculate beta in radians
        beta = -1*betaRange + j*betastep;
        % Y and Z components of coil orientation
        Wy = cos(alpha)*sin(beta);
        Wz = cos(alpha)*cos(beta);
        
        % Calculate vector from origin to resonator,
        % and then vector from resonator to receiver
        resLoc = radius.*[Wx,Wy,Wz];
        R = recLoc - resLoc;
        
        % Make R a unit vector
        Rtot = sqrt(R(1)^2 + R(2)^2 + R(3)^2);
        R = R ./ Rtot;
        
        % Calculate B-field vector based on coil orientation
        % and vector from resonator to receiver
        d = dot([Wx, Wy, Wz], R);
        B = (3/2)*d*R - ([Wx, Wy, Wz] ./ 2);
        Btot = sqrt(B(1)^2 + B(2)^2 + B(3)^2);
        % Normalize B-field components
        X = B(1)/Btot;
        Y = B(2)/Btot;
        Z = B(3)/Btot;
        
        % Make sure Z is positive
        if (Z < 0)
            X = -1*X;
            Y = -1*Y;
            Z = -1*Z;
        end
        
        % Store X component in the test matrix
        Bx((i-1)*b+j) = X;
    end
end

% Sort X component array
sortX = sort(Bx);
xdivisionsize = size/xnumdivisions;

% Pick value out of sorted X's at each multiple of division size
xdivisions = zeros(xnumdivisions, 1);
xdivisions(1) = -1;                     % Manually create lowest division
for xdiv = 2:xnumdivisions
    xdivisions(xdiv) = sortX((xdiv-1)*xdivisionsize);
end

% Empty array to store Y values
ByAdjust = zeros(xnumdivisions, xdivisionsize);

% Analysis parameter for how many Y values not placed
notplaced = 0;

% Y component distribution-finding loop
for i = 1:a
    % Calculate alpha in radians
    alpha = resAlpha - alphaRange +  i*alphastep;
    % X component of coil orientation
    Wx = sin(alpha);
    
    for j = 1:b
        % Calculate beta in radians
        beta = -1*betaRange + j*betastep;
        % Y and Z components of coil orientation
        Wy = cos(alpha)*sin(beta);
        Wz = cos(alpha)*cos(beta);
        
        % Calculate vector from origin to resonator,
        % and then vector from resonator to receiver
        resLoc = radius.*[Wx,Wy,Wz];
        R = recLoc - resLoc;
        
        % Make R a unit vector
        Rtot = sqrt(R(1)^2 + R(2)^2 + R(3)^2);
        R = R ./ Rtot;
        
        % Calculate B-field vector based on coil orientation
        % and vector from res to RX
        d = dot([Wx, Wy, Wz], R);
        B = (3/2)*d*R - ([Wx, Wy, Wz] ./ 2);
        Btot = sqrt(B(1)^2 + B(2)^2 + B(3)^2);
        % Normalize B-field components
        X = B(1)/Btot;
        Y = B(2)/Btot;
        Z = B(3)/Btot;
        
        % Make sure Z is positive
        if (Z < 0)
            X = -1*X;
            Y = -1*Y;
            Z = -1*Z;
        end
        
        index = 1;
        % Calculate the corresponding X value -> determine which
        % X division we are in
        for d = 2:xnumdivisions
            if (X > xdivisions(d)) index = index + 1;
            end
        end
        
        % Place new Y value into the adjusted Y array,
        % in the row designated by the X division,
        % in the next available spot
        currentRow = ByAdjust(index, :);
        newRow = placeNewVal(Y, currentRow);
        if (currentRow(size/xnumdivisions) ~= 0 &&...
                currentRow(size/xnumdivisions) == newRow(size/xnumdivisions))
            notplaced = notplaced + 1;
        end
        ByAdjust(index, :) = newRow;
    end
end

notplaced

% Need to fill the not-placed -1s
for xdiv = 1:xnumdivisions
    row = ByAdjust(xdiv, :);
    recurseIndex = xdivisionsize;
    ByAdjust(xdiv, :) = recurseFill(row, recurseIndex);
end

% Sort Y component array
ydivisionsize = size/xnumdivisions/ynumdivisions;
ydivisions = zeros(xnumdivisions, ynumdivisions);
sortY = ByAdjust;
for xdiv = 1:xnumdivisions
    currentRow = sort(ByAdjust(xdiv, :)); 
    sortY(xdiv, :) = currentRow;
    % Compile row of Y divisions for each X division
    newDivisionSet = zeros(1, ynumdivisions);
    newDivisionSet(1) = -1;
    for ydiv = 2:ynumdivisions
        newDivisionSet(ydiv) = currentRow((ydiv-1)*ydivisionsize);
    end
    ydivisions(xdiv, :) = newDivisionSet;
end

% Plots of Thresholds from -1 to 1
close all
zlabels = linspace(1,xnumdivisions,xnumdivisions);
xlabels = linspace(1,ynumdivisions,ynumdivisions);
scatterSize = 9;
scatter(zlabels,xdivisions, scatterSize);
title('Z Component Thresholds');
xlabel('Threshold Number (1-64)');
ylabel('Z Value');
axis([1 64 -1 1]);

figure
for ydiv = 1:xnumdivisions
    scatter(xlabels, ydivisions(ydiv, :), scatterSize);
    hold on;
end
title('X Component Thresholds by Z Division');
xlabel('Threshold Number (1-32)');
ylabel('X Value');
axis([1 32 -1 1]);

%% Loop to Fill Table

tableIndices = zeros(size, 2);

for i = 1:a
    % Calculate eyeball alpha in radians
    eyeAlpha = -1*alphaRange +  i*alphastep;
    alpha = eyeAlpha + resAlpha;
    % X component of coil orientation
    Wx = sin(alpha);
    
    for j = 1:b
        % Calculate beta in radians
        beta = -1*betaRange + j*betastep;
        % Y and Z components of coil orientation
        Wy = cos(alpha)*sin(beta);
        Wz = cos(alpha)*cos(beta);
        
        % Calculate vector from origin to resonator,
        % and then vector from resonator to receiver
        resLoc = radius.*[Wx,Wy,Wz];
        R = recLoc - resLoc;
        
        % Make R a unit vector
        Rtot = sqrt(R(1)^2 + R(2)^2 + R(3)^2);
        R = R ./ Rtot;
        
        % Calculate B-field vector based on coil orientation
        % and vector from res to RX
        d = dot([Wx, Wy, Wz], R);
        B = (3/2)*d*R - ([Wx, Wy, Wz] ./ 2);
        Btot = sqrt(B(1)^2 + B(2)^2 + B(3)^2);
        % Normalize B-field components
        X = B(1)/Btot;
        Y = B(2)/Btot;
        Z = B(3)/Btot;
        
        % Make sure Z is positive
        if (Z < 0)
            X = -1*X;
            Y = -1*Y;
            Z = -1*Z;
        end
        
        xindex = 1;
        % Calculate X contribution to index using divisions
        % determined by distribution finder above
        for d = 2:xnumdivisions
            if (X > xdivisions(d))
                xindex = xindex + 1;
            end
        end
        
        ydivisionsForXDiv = ydivisions(xindex, :);
        
        yindex = 1;
        % Calculate Y contribution to index using divisions
        % determined by distribution finder above
        for d = 2:ynumdivisions
            if (Y > ydivisionsForXDiv(d))
                yindex = yindex + 1;
            end
        end
        
        linearIndex = (i-1)*b + j;
        
        tableIndices(linearIndex, 1) = xindex;
        tableIndices(linearIndex, 2) = yindex;
        
        % Place angles
        lookup(xindex, yindex, 1, :) =...
            placeNewVal(eyeAlpha, lookup(xindex, yindex, 1, :));
        lookup(xindex, yindex, 2, :) =...
            placeNewVal(beta, lookup(xindex, yindex, 2, :));
    end
end

%% Lookup Table Export Preparation

% New lookup table array to fill
angleLookup = zeros(xnumdivisions, ynumdivisions, 2);

% Calculate average alpha and beta values for each row in original
% lookup table, and insert those into the new table
% First, alpha and beta are shifted by pi/2 to
% change range from -range:range to 0:2*range
% Alpha values are multiplied by 65,535/(2*range)
% Beta values are multiplied by 65,535/(2*range)
% Purpose of this is to convert to unsigned int
for xindex = 1:xnumdivisions
    for yindex = 1:ynumdivisions
        alphaList = lookup(xindex, yindex, 1, :);
        betaList =lookup(xindex, yindex, 2, :);
        
        alphaVal = meanNoZeros(alphaList);
        alphaVal = (alphaVal + alphaRange)*65535/(2*alphaRange);
        betaVal = meanNoZeros(betaList);
        betaVal = (betaVal + betaRange)*65535/(2*betaRange);
        
        angleLookup(xindex, yindex, 1) = alphaVal;
        angleLookup(xindex, yindex, 2) = betaVal;
    end
end

%% Export

% Export Data Types:
% X Thresholds - Signed 16-bit int
% Y Thresholds - Signed 16-bit int
% Angles - Unsigned 16-bit int

xdivisionsNew = 32767 .* xdivisions;
ydivisionsNew = 32767 .* ydivisions;

% Export X thresholds
fid = fopen('xThresh.h', 'wt');
fprintf(fid, 'const PROGMEM int16_t xThresh[%d]={\n%5.0f',...
    xnumdivisions, xdivisionsNew(1));
fprintf(fid, ',\n%5.0f',xdivisionsNew(2:end));
fprintf(fid,'};');
fclose(fid);

% Export Y thresholds
fid = fopen('yThresh.h', 'wt');
fprintf(fid, 'const PROGMEM int16_t yThresh[%d][%d]={\n%5.0f',...
    xnumdivisions, ynumdivisions, ydivisionsNew(1, 1));
fprintf(fid, ',%5.0f',ydivisionsNew(1, 2:end));
for row = 2:xnumdivisions
    fprintf(fid, ',\n%5.0f', ydivisionsNew(row, 1));
    fprintf(fid, ',%5.0f', ydivisionsNew(row, 2:end));
end
fprintf(fid,'};');
fclose(fid);

% Export lookup table
fid = fopen('angleLookup.h', 'wt');
fprintf(fid, 'const PROGMEM uint16_t angleLookup[%d][%d][%d]={\n',...
    xnumdivisions, ynumdivisions, 2);
for row = 1:(xnumdivisions-1)
    fprintf(fid, '{ ');
    for column = 1:(ynumdivisions-1)
        fprintf(fid, '{%5.0f, %5.0f}, ',...
            angleLookup(row, column, 1),...
            angleLookup(row, column, 2));
    end
    fprintf(fid, '{%5.0f, %5.0f} },\n',...
        angleLookup(row, ynumdivisions, 1),...
        angleLookup(row, ynumdivisions, 2));
end
fprintf(fid, '{ ');
for column = 1:(ynumdivisions-1)
    fprintf(fid, '{%5.0f, %5.0f}, ',...
            angleLookup(xnumdivisions, column, 1),...
            angleLookup(xnumdivisions, column, 2));
end
fprintf(fid, '{%5.0f, %5.0f} }\n};',...
    angleLookup(xnumdivisions, ynumdivisions, 1),...
    angleLookup(xnumdivisions, ynumdivisions, 2));
fclose(fid);


%% Analysis
% Determine how much of the table is filled
% and how the angles in each bucket compare to each other
num = 0;
maxlength = 0;
filledSizes = zeros(p);

% Vectors to store average alpha/beta and differences
% between max and min alpha/beta
diffVect = zeros(p, 1);
alphaVect = zeros(p, 1);
betaVect = zeros(p, 1);

for row = 1:xnumdivisions
    for column = 1:ynumdivisions
        % Extract angle vectors at row, column
        alphaList = lookup(row, column, 1, :);
        betaList = lookup(row, column, 2, :);
        
        numfilled = 0;
        
        % For each angle pair stored in bucket
        for q = 1:listlength
            alphaVal = alphaList(q);
            betaVal = betaList(q);
            
            % If there's a value, analyze it
            if (alphaVal ~= 0)
                % Add to number of filled spots in bucket
                numfilled = numfilled + 1;
                
                
            % If it's just zero, we're done
            else break;
            end
        end
        
        % Calculate 1D index from 2D indices
        linearIndex = (row-1)*ynumdivisions + column;
        
        % Place difference in vector
        diff = vectorDifference(alphaList, betaList);
        diffVect(linearIndex) = diff*180/pi;
        
        % Place average in vector
        alphaVect(linearIndex) = sum(alphaList)/numfilled*180/pi;
        betaVect(linearIndex) = sum(betaList)/numfilled*180/pi;
        
        % Update maximum bucket length
        if (numfilled > maxlength) maxlength = numfilled;
        end
        
        % Update number of buckets wtih something in them
        if (numfilled > 0) num = num + 1;
        end
        
        % Number of filled buckets placed in vector
        filledSizes(linearIndex) = numfilled;
    end
end

% Fraction of filled buckets
frac = num / p

% Max bucket length
maxlength

% Maximum orientation difference in same bucket
maxDifferenceWithinBucket = max(diffVect)

% Heatmap loop
heatX = zeros(size, 1);
heatY = zeros(size, 1);
heatZ = zeros(size, 1);
heatColor = zeros(size, 1);

for i = 1:a
    
    % Calculate eyeball alpha in radians
    eyeAlpha = -1*alphaRange +  i*alphastep;
    alpha = eyeAlpha + resAlpha;
    % X component of coil orientation
    Wx = sin(alpha);
    
    for j = 1:b
        % Calculate beta in radians
        beta = -1*betaRange + j*betastep;
        % Y and Z components of coil orientation
        Wy = cos(alpha)*sin(beta);
        Wz = cos(alpha)*cos(beta);
        
        % Calculate vector from origin to resonator,
        % and then vector from resonator to receiver
        resLoc = radius.*[Wx,Wy,Wz];
        R = recLoc - resLoc;
        
        % Make R a unit vector
        Rtot = sqrt(R(1)^2 + R(2)^2 + R(3)^2);
        R = R ./ Rtot;
        
        % Calculate B-field vector based on coil orientation
        % and vector from res to RX
        d = dot([Wx, Wy, Wz], R);
        B = (3/2)*d*R - ([Wx, Wy, Wz] ./ 2);
        Btot = sqrt(B(1)^2 + B(2)^2 + B(3)^2);
        % Normalize B-field components
        X = B(1)/Btot;
        Y = B(2)/Btot;
        Z = B(3)/Btot;
        
        % Make sure Z is positive
        if (Z < 0)
            X = -1*X;
            Y = -1*Y;
            Z = -1*Z;
        end
        
        xindex = 1;
        % Calculate X contribution to index using divisions
        % determined by distribution finder above
        for d = 2:xnumdivisions
            if (X > xdivisions(d))
                xindex = xindex + 1;
            end
        end
        
        ydivisionsForXDiv = ydivisions(xindex, :);
        
        yindex = 1;
        % Calculate Y contribution to index using divisions
        % determined by distribution finder above
        for d = 2:ynumdivisions
            if (Y > ydivisionsForXDiv(d))
                yindex = yindex + 1;
            end
        end
        
        alphaList = lookup(xindex, yindex, 1, :);
        betaList = lookup(xindex, yindex, 2, :);
        
        alphaFromLookup = meanNoZeros(alphaList);
        betaFromLookup = meanNoZeros(betaList);
        
        orientation = [sin(eyeAlpha), cos(eyeAlpha)*sin(beta), ...
            cos(eyeAlpha)*cos(beta)];
        orientationFromLookup = [sin(alphaFromLookup),...
            cos(alphaFromLookup)*sin(betaFromLookup),...
            cos(alphaFromLookup)*cos(betaFromLookup)];
        difference = acos(dot(orientation, orientationFromLookup));
        
        linearIndex = (i-1)*b + j;
        
        heatX(linearIndex) = orientation(1);
        heatY(linearIndex) = orientation(2);
        heatZ(linearIndex) = orientation(3);
        heatColor(linearIndex) = difference*180/pi;
    end
end

%% Plots

maxDifferenceActualToRepresentative = max(heatColor)

scatterSize = 9;

% Plots of difference in eye orientation values as a function of average
figure
scatter3(alphaVect, betaVect, diffVect, scatterSize);
title('Range in Orientation Sorted into each Bin');
xlabel('Average \phi (\circ)');
ylabel('Average \theta (\circ)');
zlabel('Range of Orientation (\circ)');

% Heatmap of error in orientation from lookup as a function of 
% actual orientation
figure
scatter3(heatX, heatY, heatZ, scatterSize, heatColor);
title('Heatmap: Error in Measured Orientation');
xlabel('Z');
ylabel('X');
zlabel('Y');
caxis([0 3]);
cb = colorbar;
title(cb, 'Error (\circ)');
[x,y,z] = sphere;
lightGrey = 0.8*[1 1 1];
surface(x,y,z, 'FaceColor', 'none', 'EdgeColor', lightGrey);
axis([-1 1 -1 1 0 1]);
grid off

% Histogram of same ^
figure
histogram(heatColor, 30);
title('Histogram of Orientation Error');
xlabel('Error between Measured and Actual Orientation (\circ)');
ylabel('Number of Occurances (out of 65536)');

% Cumulative distribution function of same ^
figure
cdfplot(heatColor);
title('Error in Measured Orientation: CDF');
xlabel('e');
ylabel('p(E < e)');

% Plots of bucket filling
figure
plot(filledSizes);
title('Size of Bins: Ordered by Index');
xlabel('Index');
ylabel('Number of Occupants');
axis([0 4500 0 2000]);
figure
plot(sort(filledSizes));
title('Size of Bins: Ordered by Size');
xlabel('Bin');
ylabel('Number of Occupants');

%% Helper Functions
    % Takes a number, and a number representing a power of 2,
    % and returns a bit indicating whether the power of two is smaller
    % than the number, and if so, returns the leftover portion
    % e.g. getMSB(18, 4) = [1, 2], getMSB(18, 5) = [0, 18]
    function [bit, leftover] = getMSB(num, power)
        twoPower = 2^power;
        if (num < twoPower)
            bit = 0;
            leftover = num;
        else
            bit = 1;
            leftover = num - twoPower;
        end
    end

    % Takes a 1D array, and a value, and places the value in the next
    % unfilled spot in the array
    function newList = placeNewVal(val, currentList)
        len = length(currentList);
        for k = 1:len
            if (currentList(k) == 0)
                currentList(k) = val;
                break;
            end
        end
        newList = currentList;
    end

    % Takes a row of y divisions, which may have empty zeros at the end,
    % and fills those zeros with -1
    function fullRow = recurseFill(row, index)
        if (index == 0 || row(index) ~= 0)
            fullRow = row;
        else
            row(index) = -1;
            fullRow = recurseFill(row, index-1);
        end                
    end

    % Takes a 1D array, which may contain empty zeros at the end,
    % and calculates the mean of the filled values
    function m = meanNoZeros(array)
        len = length(array);
        s = 0;
        for k = 1:len
            val = array(k);
            if (val == 0)
                m = s/(k-1);
                break;
            else
                s = s + val;
            end
        end
    end

    % Takes two lists representing alpha and beta values stored
    % in a single lookup table bucket, and computes the orientation
    % signified by each pair, returning the maximum anglular
    % difference between two such orientations in the bucket
    function ang = vectorDifference(alphaList, betaList)
        
        len = length(alphaList);
        % If there is only one value in the bucket,
        % difference is zero
        if (len == 1)
            ang = 0;
        else
            % Compute orientation vector from alpha and beta
            thisAlpha = alphaList(1);
            thisBeta = betaList(1);
            thisVector = [sin(thisAlpha),...
                cos(thisAlpha)*sin(thisBeta),...
                cos(thisAlpha)*cos(thisBeta)];
            maxAng = 0;
            for ind = 2:len
                otherAlpha = alphaList(ind);
                otherBeta = betaList(ind);
                
                % If the other angles are in fact empty
                % (0 represents unfilled slot), we're done
                if (otherAlpha == 0)
                    angleBetween = 0;
                else
                    % Otherwise, compute orientation vector
                    otherVector = [sin(otherAlpha),...
                        cos(otherAlpha)*sin(otherBeta),...
                        cos(otherAlpha)*cos(otherBeta)];
                    % Use dot product to find difference
                    angleBetween = acos(dot(thisVector, otherVector));
                end
                % If difference exceeds previous max, it is new max
                if (angleBetween > maxAng)
                    maxAng = angleBetween;
                end
            end
            
            % Recursively check every other angle pair in list
            angFromRest = vectorDifference(alphaList(2:len), betaList(2:len));
            if (angFromRest > maxAng)
                ang = angFromRest;
            else
                ang = maxAng;
            end
        end
    end

end