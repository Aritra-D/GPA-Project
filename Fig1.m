% Fig 1

function Fig1(combineUniqueElectrodeData,getSpikeElectrodesFlag,unitID,spikeCutoff,snrCutoff,ECoGFlag)


if ~exist('combineUniqueElectrodeData','var'); combineUniqueElectrodeData=1;     end
if ~exist('getSpikeElectrodesFlag','var');     getSpikeElectrodesFlag=1;         end  % Set to 1 if you want to analyze spike electrodes
if ~exist('unitID','var');                     unitID=1;                         end  % 0 - unsorted hash. 1 - sorted Spikes
if ~exist('spikeCutoff','var');                spikeCutoff=25;                   end
if ~exist('snrCutoff','var');                  snrCutoff=2;                      end
if ~exist('ECoGFlag','var');                   ECoGFlag=0;                       end % 0 = all, 1=ECoG only, 2==microECoG only;

%%%%%%%%%%%%%%%%%%%%%%%%%%% display properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%
fontSizeSmall=10; fontSizeMedium=12; fontSizeLarge=14;
backgroundColor = 'w'; labelSize=10;
attOutColors   = repmat('k',8,1);  attInColors    = jet(8);

MonkeyNameLists{1} = 'abu'; MonkeyNameLists{2} = 'rafiki'; gridType = 'Microelectrode';
freqRanges{1} = [8 12];
freqRanges{2} = [30 58];
freqRanges{3} = [102 238];
freqRanges{4} = [250 500];
freqRangeStr = {'alpha','gamma','hi-gamma','sigma'};
numFreqRanges = length(freqRanges);

for MonkeyIndex = 1: size(MonkeyNameLists,2)
    clear monkeyName;
    monkeyName = MonkeyNameLists{MonkeyIndex};
%%%%%%%%%%%%%%%%%%%%%%%%%% Find Good Electrodes %%%%%%%%%%%%%%%%%%%%%%%%%%%
[allGoodElectrodes,allDaysToCombine,~] = getGoodElectrodes(monkeyName,gridType,combineUniqueElectrodeData,getSpikeElectrodesFlag,unitID,spikeCutoff,snrCutoff,ECoGFlag);
[expDates,protocolNames,~,folderSourceString] = dataInformation(MonkeyNameLists{MonkeyIndex},gridType,ECoGFlag);
numGoodElectrodes = length(allGoodElectrodes);

%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Handles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(monkeyName, 'abu')
hDEnergyVsFreq       = getPlotHandles(1,3,[0.15 0.8 0.7 0.18],0.03,0.01);
hBarPlot            = getPlotHandles(1,3,[0.15 0.55 0.7 0.18],0.03,0.01);

elseif strcmpi(monkeyName, 'rafiki')
hDEnergyVsFreq       = getPlotHandles(1,3,[0.15 0.30 0.7 0.18],0.03,0.01);
hBarPlot            = getPlotHandles(1,3,[0.15 0.05 0.7 0.18],0.03,0.01);
end

% LFP Power is always computed for this period. ERP and Firing Rates are
% computed for the range timeRangeForComputation.

timeRangeForBLComputation = [-0.25 -0.05];
timeRangeForSTComputation = [0.2 0.4];

    clear energyDataTMP energyDataBLTMP 
    for i=1:numGoodElectrodes
        disp([num2str(i) ' of ' num2str(numGoodElectrodes)]);
        [energyDataTMP(i),energyDataBLTMP(i)] = getEnergyDataMT(monkeyName,expDates,protocolNames,folderSourceString,gridType,...
                    allGoodElectrodes(i),allDaysToCombine{i},timeRangeForSTComputation,timeRangeForBLComputation,freqRanges);
        for j=1:numFreqRanges
            meanE{j}(i,:,:) = energyDataTMP(i).analysisData{j};
        end
                
            for aPos=1:2
                for cPos=6:8
                    dEnergyVsFrequency(i,aPos,cPos,:) = squeeze(energyDataTMP(i).data(aPos,cPos,:))' - squeeze(energyDataBLTMP(i).data(1,cPos,:))'; % Note 1 instead of aPos for BL
                    dEnergyVsFrequencyBL(i,aPos,cPos,:) = squeeze(energyDataBLTMP(i).data(aPos,cPos,:))' - squeeze(energyDataBLTMP(i).data(1,cPos,:))'; % Note 1 instead of aPos for BL
                end
            end
    end
            % combine Data
            energyData = combineEnergyData(energyDataTMP);
            energyDataBL = combineEnergyData(energyDataBLTMP);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for cPos=6:8
            %%%%%%%%%%%%%%%%%%%%%%% EnergyVsFreq %%%%%%%%%%%%%%%%%%%%%%%%%%
%             hold(hEnergyVsFreq(1,cPos),'on');
%             plot(hEnergyVsFreq(1,cPos),energyDataBL.freqVals,squeeze(energyDataBL.data(1,cPos,:)),'color', [0.8 0.8 0.8]);
%             plot(hEnergyVsFreq(1,cPos),energyDataBL.freqVals,squeeze(energyDataBL.data(2,cPos,:)),'color', [0.8 0.8 0.8]);
%             
%             plot(hEnergyVsFreq(1,cPos),energyData.freqVals,squeeze(energyData.data(1,cPos,:)),'color', attOutColors(cPos,:));
%             plot(hEnergyVsFreq(1,cPos),energyData.freqVals,squeeze(energyData.data(2,cPos,:)),'color', attInColors(cPos,:));
            figure(1)
%             
            %%%%%%%%%%%%%%%%%%%%%% DEnergyVsFreq %%%%%%%%%%%%%%%%%%%%%%%%%%
            hold(hDEnergyVsFreq(1,cPos-5),'on');
            %plot(hDEnergyVsFreq(1,cPos),energyData.freqVals,zeros(1,length(energyData.freqVals)),'color',[0.8 0.8 0.8]);
            plot(hDEnergyVsFreq(1,cPos-5),energyData.freqVals,squeeze(mean(dEnergyVsFrequency(:,1,cPos,:),1)),'color', attOutColors(cPos,:),'LineWidth',2);
            plot(hDEnergyVsFreq(1,cPos-5),energyData.freqVals,squeeze(mean(dEnergyVsFrequency(:,2,cPos,:),1)),'color', attInColors(cPos,:),'LineWidth',2);
            plot(hDEnergyVsFreq(1,cPos-5),energyData.freqVals,squeeze(mean(dEnergyVsFrequencyBL(:,1,cPos,:),1)),'color', [0.8 0.8 0.8]);
            plot(hDEnergyVsFreq(1,cPos-5),energyData.freqVals,squeeze(mean(dEnergyVsFrequencyBL(:,2,cPos,:),1)),'color', [0.4 0.4 0.4]);
            
            if cPos==6
               xlabel(hDEnergyVsFreq(1,cPos-5),'Frequency(Hz)');
               ylabel(hDEnergyVsFreq(1,cPos-5),'log(Power)');
%             else   
%                set(hDEnergyVsFreq(1,cPos-5),'YTickLabel',[],'fontSize',labelSize);
            end
            
            
            
            displayMeanBarPlot(hBarPlot(1,1),meanE{1});
            displayMeanBarPlot(hBarPlot(1,2),meanE{2});
            
            % Peak and Max Gamma
            [peakGamma,maxGamma] = getPeakFrequency(dEnergyVsFrequency,energyData.freqVals);
            displayMeanBarPlot(hBarPlot(1,3),peakGamma);
            
            if cPos==8
                xlabel(hBarPlot(1,cPos-5),'Contrast(%)');
                ylabel(hBarPlot(1,cPos-5),'Gamma Peak Frequncy');
            else
                xlabel(hBarPlot(1,cPos-5),'Contrast(%)');
                ylabel(hBarPlot(1,cPos-5),'log(Power)');
            end
            
            figure(1)
            end
            
            
            
%             yticklabels(hBarPlot(1,3),{ '0', '20', '40','60'});
            
end
end

function [energyData,energyDataBL] = getEnergyDataMT(monkeyName,expDates,protocolNames,folderSourceString,gridType,electrodeList,dayList,timeRangeForSTComputation,timeRangeForBLComputation,freqRanges)

numElectrodes = length(dayList);

    if length(electrodeList)==1
        electrodeList = repmat(electrodeList,[1 numElectrodes]);
    end
    
    clear erpDataTMP firingRateDataTMP energyDataTMP energyDataBLTMP
    for i=1:numElectrodes
        [energyDataTMP(i),energyDataBLTMP(i)] = getEnergyDataMTSingleDay(monkeyName,expDates{dayList(i)},protocolNames{dayList(i)},folderSourceString,gridType,electrodeList(i),timeRangeForSTComputation,timeRangeForBLComputation,freqRanges);
        
          % No need to normalize energy because we use log scale
        energyData   = combineEnergyData(energyDataTMP);
        energyDataBL = combineEnergyData(energyDataBLTMP);
    
    end
end
function [energyData,energyDataBL] = getEnergyDataMTSingleDay(monkeyName,expDate,protocolName,folderSourceString,gridType,electrodeNum,timeRangeForSTComputation,timeRangeForBLComputation,freqRanges)



[parameterCombinations,cValsUnique] = loadParameterCombinations(monkeyName,expDate,protocolName,folderSourceString,gridType); 
numContrasts = length(cValsUnique);

% Get bad trials
load([folderSourceString 'data\' monkeyName '\' gridType '\' expDate '\' protocolName '\segmentedData\badTrials.mat']);

% Get timeVals
load([folderSourceString 'data\' monkeyName '\' gridType  '\' expDate '\' protocolName '\segmentedData\LFP\lfpInfo.mat']);


% Get positions for MT computations
Fs              =  round(1/(timeVals(2)-timeVals(1)));
params.tapers   = [2 3];
params.pad      = -1;
params.Fs       = Fs;
params.fpass    = [0 100];
params.trialave = 1;
    
stPos  = find(timeVals>=timeRangeForSTComputation(1),1) + (0:round(Fs*diff(timeRangeForSTComputation))-1);
blPos  = find(timeVals>=timeRangeForBLComputation(1),1) + (0:round(Fs*diff(timeRangeForBLComputation))-1);

% Get Analog Data
load([folderSourceString 'data\' monkeyName '\' gridType  '\' expDate '\' protocolName '\segmentedData\LFP\elec' num2str(electrodeNum) '.mat']);
        
N = zeros(2,numContrasts);
numFreqs = length(freqRanges);

for aPos=1:2 % 1-out, 2-in
    for cPos=6:numContrasts
        clear goodPos
        goodPos = setdiff(parameterCombinations{cPos,1,1,aPos,1},badTrials);
        N(aPos,cPos) = length(goodPos);
            
             % Multi-tapering
                % Power
                data       = analogData(goodPos,stPos)';
                tmpEST     = mtspectrumc(data,params);
                data       = analogData(goodPos,blPos)';
                [tmpEBL,ys]= mtspectrumc(data,params);
                
                [freqPosToRemove,frequenciesToRemove] = getFreqPosToRemove(ys); %#ok<NASGU>
                
                mEnergyVsFreqST(aPos,cPos,:) = conv2Log(tmpEST);
                mEnergyVsFreqBL(aPos,cPos,:) = conv2Log(tmpEBL);
                %mEnergyVsFreqST(aPos,cPos,:) = conv2Log(removeNoiseAndSmooth(tmpEST(:),ys,frequenciesToRemove));
                %mEnergyVsFreqBL(aPos,cPos,:) = conv2Log(removeNoiseAndSmooth(tmpEBL(:),ys,frequenciesToRemove));
                
                for i=1:numFreqs
                    energyVals{i}(aPos,cPos) = conv2Log(getMeanEnergyForAnalysis(tmpEST(:),ys,freqRanges{i},freqPosToRemove)); %
                    energyValsBL{i}(aPos,cPos) = conv2Log(getMeanEnergyForAnalysis(tmpEBL(:),ys,freqRanges{i},freqPosToRemove)); %,
                end
    end
end



    % Multi-tapering
    energyData.data = mEnergyVsFreqST;
    energyData.analysisData = energyVals;
    energyData.freqVals = ys;
    energyData.N = N;
    
    energyDataBL.data = mEnergyVsFreqBL;
    energyDataBL.analysisData = energyValsBL;
    energyDataBL.freqVals = ys;
    energyDataBL.N = N;


end
function combinedX = combineData(x)

numEntries = length(x);

% Initialize
combinedX.data = zeros(size(x(1).data));
combinedX.analysisData = zeros(size(x(1).analysisData));
combinedX.timeVals = x(1).timeVals;
combinedX.N = 0;

for i=1:numEntries
    combinedX.data = combinedX.data + (x(i).data)/numEntries;
    combinedX.analysisData = combinedX.analysisData + (x(i).analysisData)/numEntries;
    combinedX.N = combinedX.N + x(i).N;
end
end
function combinedX = combineEnergyData(x)

numEntries = length(x);

% Initialize
combinedX.data = zeros(size(x(1).data));
numFreqRanges = size(x(1).analysisData,2);
combinedX.analysisData = cell(1,numFreqRanges);
for j=1:numFreqRanges
    combinedX.analysisData{j} = zeros(size(x(1).analysisData{1}));
end
combinedX.freqVals = x(1).freqVals;
combinedX.N = 0;

for i=1:numEntries
    combinedX.data = combinedX.data + (x(i).data)/numEntries;
    for j=1:numFreqRanges
        combinedX.analysisData{j} = combinedX.analysisData{j} + (x(i).analysisData{j})/numEntries;
    end
    combinedX.N = combinedX.N + x(i).N;
end
end
function mEnergyOut = removeNoiseAndSmooth(mEnergy,freq,frequenciesToRemove)

if ~isempty(frequenciesToRemove)
    N = size(frequenciesToRemove,1);
    for i=1:N
        fPosMin = getClosestIndex(freq,frequenciesToRemove(i,1)); fPosMax = getClosestIndex(freq,frequenciesToRemove(i,2));
        mEnergyMin = mEnergy(fPosMin-1); mEnergyMax = mEnergy(fPosMax+1);
        
        for j=fPosMin:fPosMax
            mEnergy(j) = mEnergyMin+ (j-fPosMin-1)*(mEnergyMax-mEnergyMin)/(fPosMax-fPosMin+2);
        end
    end
end

mEnergyOut = smooth(mEnergy);
end
function index = getClosestIndex(xs,val)
xs2 = abs(xs-val);
index = find(xs2==min(xs2),1);
end
function [freqPosToRemove,frequenciesToRemove] = getFreqPosToRemove(freq)

frequenciesToRemove = [58 62; 98 102; 118 122; 238 242];

N = size(frequenciesToRemove,1);
freqPosToRemove=[];
for i=1:N
    freqPosToRemove = cat(2,freqPosToRemove,getClosestIndex(freq,frequenciesToRemove(i,1)):getClosestIndex(freq,frequenciesToRemove(i,2)));
end
end

function [peakGamma,maxGamma] = getPeakFrequency(dEnergyVsFrequency,ys)

gammaRange = [20 60];
gammaPos = intersect(find(ys>=gammaRange(1)),find(ys<gammaRange(2)));

numEntries = size(dEnergyVsFrequency,1);

for i=1:numEntries
    for aPos=1:2
        for cPos=1:8
            x = squeeze(dEnergyVsFrequency(i,aPos,cPos,gammaPos));
            maxGamma(i,aPos,cPos) = max(x);
            peakGamma(i,aPos,cPos) = ys(gammaPos(find(x==maxGamma(i,aPos,cPos),1)));
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load functions
function [parameterCombinations,cValsUnique,tValsUnique,eValsUnique,aValsUnique,sValsUnique] = loadParameterCombinations(monkeyName,expDate,protocolName,folderSourceString,gridType) %#ok<STOUT>
load([folderSourceString 'data\' monkeyName '\' gridType '\' expDate '\' protocolName '\extractedData\parameterCombinations.mat']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function eValue = getMeanEnergyForAnalysis(mEnergy,freq,freqRange,freqPosToRemove)

posToAverage = setdiff(intersect(find(freq>=freqRange(1)),find(freq<freqRange(2))),freqPosToRemove);
eValue   = mean(mEnergy(posToAverage));
end
function displayMeanBarPlot(hPlot,meanData)

numContrasts=8;
contrastIndices = 0:numContrasts-1;

meanDataOut = squeeze(meanData(:,1,:));
meanDataIn  = squeeze(meanData(:,2,:));

% Plot raw Data
numDays = size(meanData,1);

if numDays>1
    mOut = squeeze(mean(meanDataOut,1)); semOut = squeeze(std(meanDataOut,0,1))/sqrt(numDays);
    mIn  = squeeze(mean(meanDataIn,1));  semIn  = squeeze(std(meanDataIn,0,1))/sqrt(numDays);
else
    mOut = meanDataOut; semOut = zeros(1,numContrasts);
    mIn  = meanDataIn;  semIn  = zeros(1,numContrasts);
end

hold(hPlot,'on');
i = 6;
AttInOutBar = [mOut(i) mIn(i); mOut(i+1) mIn(i+1); mOut(i+2) mIn(i+2);];
semAttInOut = [semOut(i) semIn(i); semOut(i+1) semIn(i+1); semOut(i+2) semIn(i+2);];
bar(hPlot,AttInOutBar)

% Properties of the bar graph as required
% ax.YGrid = 'on';
% ax.GridLineStyle = '-';
xticks(hPlot,[1 2 3]);

% Naming each of the bar groups
xticklabels(hPlot,{ '25%', '50%', '100%'});

% X and Y labels
% xlabel (hPlot,('Contrast'));
% ylabel (hPlot,('log(Power)'));

% % Creating a legend and placing it outside the bar plot
% lg = legend('Attend-Out','Attend-In','AutoUpdate','off');
% lg.Location = 'BestOutside';
% lg.Orientation = 'Horizontal';

hold on;

% Finding the number of groups and the number of bars in each group
ngroups = size(AttInOutBar, 1);
nbars = size(AttInOutBar, 2);

% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));

% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for j = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*j-1) * groupwidth / (2*nbars);
    errorbar(hPlot,x, AttInOutBar(:,j), semAttInOut(:,j), 'k', 'linestyle', 'none');
end

% for i=6:numContrasts
%     bar(hPlot,mOut(i),semOut(i),'color',attOutColors(i,:),
%     
% end

end