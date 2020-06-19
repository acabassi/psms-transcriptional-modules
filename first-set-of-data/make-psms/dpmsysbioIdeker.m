
%% Generate PSM of Harbison data

cd '~/OneDrive - University of Cambridge, MRC Biostatistics Unit/PHD-PROJECTS/psms-transcriptional-modules/first-set-of-data/make-psms/'
addpath '~/OneDrive - University of Cambridge, MRC Biostatistics Unit/PHD-PROJECTS/DPMSysBio'
addpath '~/OneDrive - University of Cambridge, MRC Biostatistics Unit/PHD-PROJECTS/psms-transcriptional-modules/first-set-of-data/data/'

%%
fileName  = 'GalactoseData.csv';

samplingFreq     = 5;
nSamples         = 20000;
dataType         = 'Multinomial';
verbose          = false;
initialise       = true;
drawFigures      = false;
uniqueIdentifier = 1;
inputSeed        = 100;
gammaPrior       = [2 4];
hyperParameterSamplingFrequency = 1;

%%
DPM(fileName, uniqueIdentifier, nSamples, dataType, ...
    drawFigures, hyperParameterSamplingFrequency, verbose, initialise,...
    gammaPrior, samplingFreq, inputSeed)