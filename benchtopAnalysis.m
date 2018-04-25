% Caleb R Tulloss
% Brown Neuromotion Lab
% Eye Orientation Tracker
% Benchtop Data Analysis
% Version 3
% 4/17/18

function benchtopAnalysis
close all
%% data collection settings

readingsPerSample = 32;
accelPerSample = 64;

%% plot settings

scatterSize = 2;
red = [0.9, 0.11, 0.11];
green = [0.05, 0.8, 0.05];
blue = [0.11, 0.11, 0.9];

%% Time-varying analysis

% Data format
% thetaCount;phiCount;accelX;accelY;accelZ;ch1;ch2;ch3;ch4;time;temp;pressure;humid

data = dlmread('2_26_2018_V2_stable2_edited.csv');

data = data(150:end,:);

[numSamples,~] = size(data);

orientations = zeros(numSamples, 3);

% time vector (in seconds)
time = data(:,10) ./ 1000;
accelX = data(:,3);
accelY = data(:,4);
accelZ = data(:,5);
ch1 = data(:,6);
ch2 = data(:,7);
ch3 = data(:,8);
% ch4 = data(:,9);
temp = data(:,11);
pressure = data(:,12);
humid = data(:,13);

figure
plot(time, 3.3 /4096 * ch1 / readingsPerSample, 'Color', blue);
hold on
plot(time, 3.3 /4096 * ch2 / readingsPerSample, 'Color', red);
plot(time, 3.3 /4096 * ch3 / readingsPerSample, 'Color', green);
%plot(time, ch4 / readingsPerSample);
title('Output Voltage vs. Time');
xlabel('Time (s)');
ylabel('Detection Output (V)');
legend('X Position, X Coil', 'Z Position, Y Coil', 'Z Position, Z Coil');

figure
plot(time, accelX / accelPerSample);
hold on
plot(time, accelY / accelPerSample);
plot(time, accelZ / accelPerSample);
title('Output of accelerometer vs. Time');
xlabel('Time (s)');
ylabel('Accelerometer output (LSB)');
legend('X', 'Y', 'Z');

%% Weather Analysis

figure
plot(time, temp);
figure
plot(time, pressure);
figure
plot(time, humid);
title('Humidity vs. Time');
xlabel('Time (s)');
ylabel('Hmidity (%)');

voltage1 = 3.3 / 4096 * ch1 / readingsPerSample;
voltage2 = 3.3 / 4096 * ch2 / readingsPerSample;
voltage3 = 3.3 / 4096 * ch3 / readingsPerSample;

pfit1 = polyfit(humid, voltage1, 1);
pfit2 = polyfit(humid, voltage2, 1);
pfit3 = polyfit(humid, voltage3, 1);

y1 = polyval(pfit1, humid);
y2 = polyval(pfit2, humid);
y3 = polyval(pfit3, humid);

SSresid1 = sum((voltage1 - y1).^2);
SSresid2 = sum((voltage2 - y2).^2);
SSresid3 = sum((voltage3 - y3).^2);

SStotal1 = (length(voltage1)-1) * var(voltage1);
SStotal2 = (length(voltage2)-1) * var(voltage2);
SStotal3 = (length(voltage3)-1) * var(voltage3);

rsq1 = 1 - SSresid1/SStotal1
rsq2 = 1 - SSresid2/SStotal2
rsq3 = 1 - SSresid3/SStotal3

figure
hold on
scatter(humid, voltage1, scatterSize, '.', 'MarkerEdgeColor', [0.11, 0.11, 0.9]);
plot(humid, y1, 'Color', 'black');
scatter(humid, voltage2, scatterSize, '.', 'MarkerEdgeColor', [0.9, 0.11, 0.11]);
plot(humid, y2, 'Color', 'black');
scatter(humid, voltage3, scatterSize, '.', 'MarkerEdgeColor', [0.05, 0.8, 0.05]);
plot(humid, y3, 'Color', 'black');
legend('X Position, X Coil', 'Linear Fit', 'Z Position, Y Coil',...
    'Linear Fit', 'Z Position, Z Coil', 'Linear Fit');
xlabel('Humidity (%)');
ylabel('Detection Output (V)');
title('Output Voltage vs. Humidity: Linear Regression');

mean1 = 3.3 / 4096 * mean(ch1 / readingsPerSample)
mean2 = 3.3 / 4096 * mean(ch2 / readingsPerSample)
mean3 = 3.3 / 4096 * mean(ch3 / readingsPerSample)
% mean4 = mean(ch4 / readingsPerSample)
std1 = 3.3 / 4096 * std(ch1 / readingsPerSample)
std2 = 3.3 / 4096 * std(ch2 / readingsPerSample)
std3 = 3.3 / 4096 * std(ch3 / readingsPerSample)
% std4 = std(ch4 / readingsPerSample)

meanX = mean(accelX / accelPerSample)
meanY = mean(accelY / accelPerSample)
meanZ = mean(accelZ / accelPerSample)
stdX = std(accelX / accelPerSample)
stdY = std(accelY / accelPerSample)
stdZ = std(accelZ / accelPerSample)

%% Heatmap analysis

heatX = zeros(numSamples, 3);
heatY = zeros(numSamples, 3);
heatZ = zeros(numSamples, 3);

phiVector = zeros(numSamples, 1);
thetaVector = zeros(numSamples, 1);

output = zeros(numSamples, 3);

maxX = 0;
maxY = 0;
maxZ = 0;
% max4 = 0;
minX = 10000;
minY = 10000;
minZ = 10000;
% min4 = 10000;

for i = 1:numSamples
    
    accelSumX = accelX(i);
    accelSumY = accelY(i);
    accelSumZ = accelZ(i);
    
    detectionSum1 = ch1(i);
    detectionSum2 = ch2(i);
    detectionSum3 = ch3(i);
    %detectionSum4 = ch4(i);
    
    accelTotal = sqrt(accelSumX^2 + accelSumY^2 + accelSumZ^2);
    unitX = accelSumX / accelTotal;
    unitY = accelSumY / accelTotal;
    unitZ = accelSumZ / accelTotal;
    
    phi = asin(unitX);
    theta = asin(-1 * unitY / cos(phi));
    
    phiVector(i) = phi;
    thetaVector(i) = theta;
    
    orX = -1 * cos(theta) * sin(phi);
    orY = cos(theta) * cos(phi);
    orZ = -1 * sin(theta);
    
    adc1 = detectionSum1 / readingsPerSample;
    adc2 = detectionSum2 / readingsPerSample;
    adc3 = detectionSum3 / readingsPerSample;
    %adc4 = detectionSum4 / readingsPerSample;
    
    detectVoltage1 = adc1 / 4096 * 3.3;
    detectVoltage2 = adc2 / 4096 * 3.3;
    detectVoltage3 = adc3 / 4096 * 3.3;
    %detectVoltage4 = adc4 / 4096 * 3.3;
    
    output(i, :) = [detectVoltage1, detectVoltage2, detectVoltage3];
    orientations(i,:) = [orX, orY, orZ];
    heatX(i, :) = [0, 0, detectVoltage1];
    heatY(i, :) = [detectVoltage2, 0, 0];
    heatZ(i, :) = [0, detectVoltage3, 0];
    %heat4(i,:) = [detectVoltage4, 0, 0];
    
    if (detectVoltage1 > maxX) maxX = detectVoltage1;
    end
    if (detectVoltage2 > maxY) maxY = detectVoltage2;
    end
    if (detectVoltage3 > maxZ) maxZ = detectVoltage3;
    end
%     if (detectVoltage4 > max4) max4 = detectVoltage4;
%     end
    if (detectVoltage1 < minX) minX = detectVoltage1;
    end
    if (detectVoltage2 < minY) minY = detectVoltage2;
    end
    if (detectVoltage3 < minZ) minZ = detectVoltage3;
    end
%     if (detectVoltage4 < min4) min4 = detectVoltage4;
%     end
end

%% heatmap analysis

heatX = (heatX - [0, 0, minX]) ./ (maxX - minX);
heatY = (heatY - [minY, 0, 0]) ./ (maxY - minY);
heatZ = (heatZ - [0, minZ, 0]) ./ (maxZ - minZ);
%heat4 = (heat4 - [min4, 0, 0]) ./ (max4 - min4);

lightGrey = 0.8*[1 1 1];

figure
scatter3(orientations(:,3), orientations(:,1), orientations(:,2),...
    scatterSize, heatX);
title('X Channel Output Voltage as a function of Orientation');
xlabel('Z');
ylabel('X');
zlabel('Y');
[x,y,z] = sphere;
surface(x,y,z, 'FaceColor', 'none', 'EdgeColor', lightGrey);
grid off

figure
scatter3(orientations(:,3), orientations(:,1), orientations(:,2),...
    scatterSize, heatY);
title('Y Channel Output Voltage as a function of Orientation');
xlabel('Z');
ylabel('X');
zlabel('Y');
[x,y,z] = sphere;
surface(x,y,z, 'FaceColor', 'none', 'EdgeColor', lightGrey);
grid off

figure
scatter3(orientations(:,3), orientations(:,1), orientations(:,2),...
    scatterSize, heatZ);
title('Z Channel Output Voltage as a function of Orientation');
xlabel('Z');
ylabel('X');
zlabel('Y');
[x,y,z] = sphere;
surface(x,y,z, 'FaceColor', 'none', 'EdgeColor', lightGrey);
grid off

% figure
% scatter3(orientations(:,3), orientations(:,1), orientations(:,2),...
%     scatterSize, heat4);
% title('4 Channel Output Voltage as a function of Orientation');
% xlabel('Z');
% ylabel('X');
% zlabel('Y');
% [x,y,z] = sphere;
% surface(x,y,z, 'FaceColor', 'none', 'EdgeColor', lightGrey);
% grid off

heatTotal = heatX + heatY + heatZ;

figure
scatter3(orientations(:,3), orientations(:,1), orientations(:,2),...
    scatterSize, heatTotal(:, :));
title('Detection Output (RGB) vs. Orientation');
xlabel('Z');
ylabel('X');
zlabel('Y');
[z,y,x] = sphere;
surface(x,y,z, 'FaceColor', 'none', 'EdgeColor', lightGrey);
grid off
%axis([-0.6 0.6 -0.6 0.6 0 1]);
% stable V1
% axis([0.01 0.02 -0.215 -0.205 0 1]);

figure
hold on
scatter(phiVector, output(:,2), '.', 'MarkerEdgeColor', red);
scatter(phiVector, output(:,3), '.', 'MarkerEdgeColor', green);
scatter(phiVector, output(:,1), '.', 'MarkerEdgeColor', blue);
title('Output Voltage vs. \phi');
xlabel('\phi (\circ)');
ylabel('Output (V)');
legend('Y coil', 'Z coil', 'X coil');
axis([-0.2 0.3 0 3.4]);

figure
hold on
scatter(thetaVector, output(:,2), '.', 'MarkerEdgeColor', red);
scatter(thetaVector, output(:,3), '.', 'MarkerEdgeColor', green);
scatter(thetaVector, output(:,1), '.', 'MarkerEdgeColor', blue);
title('Output Voltage vs. \theta');
xlabel('\theta (\circ)');
ylabel('Output (V)');
legend('Y coil', 'Z coil', 'X coil');
axis([-0.4 0.4 0 3.4]);

figure
scatter3(thetaVector, phiVector, heatX(:,3));
figure
scatter3(thetaVector, phiVector, heatY(:,1));
figure
scatter3(thetaVector, phiVector, heatZ(:,2));

%% similarity analysis

downsampleFactor = 100;

orientations = downsample(orientations, downsampleFactor);
output = downsample(output, downsampleFactor);

numSamples = length(orientations);

numSimilarities = numSamples * (numSamples - 1) / 2

orientationAngleVSdetectionDist = zeros(numSimilarities, 5);

linearIndex = 1;

for i = 1:numSamples-1
    for j = i+1:numSamples
        orI = orientations(i, :);
        orJ = orientations(j, :);
        orAngle = radtodeg(acos(dot(orI, orJ)));
        detectI = output(i, :);
        detectJ = output(j, :);
        detectDistX = abs(detectI(1) - detectJ(1));
        detectDistY = abs(detectI(2) - detectJ(2));
        detectDistZ = abs(detectI(3) - detectJ(3));
        detectDist = sqrt(detectDistX^2 + detectDistY^2 + detectDistZ^2);
        orientationAngleVSdetectionDist(linearIndex, :) =...
            [orAngle, detectDistX, detectDistY, detectDistZ, detectDist];
        linearIndex = linearIndex + 1;
    end
end

figure
scatter(orientationAngleVSdetectionDist(:,1),...
    orientationAngleVSdetectionDist(:,5), scatterSize);
title('Difference in Detection Output vs. Difference in Orientation');
xlabel('Orientation Difference (\circ)');
ylabel('Output Difference (V)');

% sort these data by angular difference in orientation
sortedOrVSdetect = sortrows(orientationAngleVSdetectionDist, 1);

figure
plot(sortedOrVSdetect(:,1), sortedOrVSdetect(:,5));
title('Difference in Detection Output vs. Difference in Orientation');
xlabel('Orientation Difference (\circ)');
ylabel('Output Difference (V)');

numBuckets = 2230;
orientationsPerBucket = numSimilarities / numBuckets;

meanAndSTforPlot = zeros(numBuckets*2, 4);

for i = 1:numBuckets
    orAngle1 = sortedOrVSdetect((i-1)*orientationsPerBucket+1, 1);
    orAngle2 = sortedOrVSdetect(i*orientationsPerBucket, 1);
    differences = sortedOrVSdetect((i-1)*orientationsPerBucket+1:...
        i*orientationsPerBucket, 5);
    averageDifferences = mean(differences);
    stDifferences = std(differences);
    meanAndSTforPlot(i*2-1, :) = [orAngle1, averageDifferences,...
        averageDifferences + stDifferences, averageDifferences - stDifferences];
    meanAndSTforPlot(i*2, :) = [orAngle2, averageDifferences,...
        averageDifferences + stDifferences, averageDifferences - stDifferences];
end

figure
hold on
plot(meanAndSTforPlot(:,1), meanAndSTforPlot(:,2));
plot(meanAndSTforPlot(:,1), meanAndSTforPlot(:,3), 'Color', lightGrey);
plot(meanAndSTforPlot(:,1), meanAndSTforPlot(:,4), 'Color', lightGrey);
xlabel('Orientation Difference (\circ)');
ylabel('Bin-Averaged Output Difference (V)');
title('Bin-Averaged Output Difference vs. Orientation Difference');
legend('\mu', '\mu + \sigma');

end