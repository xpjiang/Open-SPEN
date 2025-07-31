%% Open SPEN in Pulseq - SPEN-SE
% Andreas Holl
% Division of Medical Physics, Department of Diagnostic and Interventional Radiology,
% University Medical Center Freiburg, Faculty of Medicine, University of Freiburg, Freiburg, Germany
% Email: andreas.holl@uniklinik-freiburg.de
% March. 23, 2024

clear all
close all
clc
%addpath(genpath([fileparts(fileparts((pwd))) '\']))

% Define FOV and resolution
fov = 300e-3;%256e-3 % 300mm
Nx = 100;
Ny = 100;
deltak = 1/fov; % Pulseq toolbox defaults to k-space units of m^-1
kWidth = Nx*deltak;

sliceThickness = 5e-3; % 5mm

% Define sequence parameters
TE = 12e-3; % 12ms
TR = 500e-3; % 20ms

rf_dur = 4e-3; % Texc=Tacq 4ms
sweepBw=25000; % 25KHz

dur_ref = 3e-3; % 3ms

% set system limits
sys = mr.opts('MaxGrad',20,'GradUnit','mT/m',...
    'MaxSlew',40,'SlewUnit','T/m/s',...
    'rfRingdownTime', 20e-6, 'rfDeadtime', 100e-6,'B0',0.3);

% Create a new sequence object
seq = mr.Sequence(sys);

% Calculate SPEN-conditions: 
% sweepBw*rf_dur = gexc.amplitude*rf_dur*fov = gacq.amplitude*Ny*mr.calcDuration(gacq)*fov



gexc = mr.makeTrapezoid('x',sys,'Amplitude',sweepBw/fov,'FlatTime',rf_dur,'Delay',sys.rfDeadTime);

% Create chirped RF and excitation gradient
rf = makeChirpedRfPulse('duration',rf_dur,'delay',sys.rfDeadTime+gexc.riseTime,'bandwidth',sweepBw, ...
    'ang',90,'n_fac',40,'system',sys); %+gexc.riseTime nachträglich geändert. Messungen waren ohne-> leichte Verzerrungen
%[bw,f0,M_xy_sta,F1]=mr.calcRfBandwidth(rf);
% gexc = mr.makeTrapezoid('x',sys,'Amplitude',bw/fov,'FlatTime',rf_dur,'Delay',sys.rfDeadTime);

% Create a slice selective sinc 180° RF refocusing pulse and crusher

[rfref, gz] = mr.makeSincPulse(pi,sys,'Duration',dur_ref,...
    'SliceThickness',sliceThickness,'apodization',0.2,'timeBwProduct',4,'PhaseOffset',0,'use','refocusing');
gzCrush = mr.makeTrapezoid('z', 'Amplitude', -gz.amplitude*3, 'Duration', (dur_ref)/2);

% Create other gradients and events
gRO = mr.makeTrapezoid('x',sys,'Flatarea',Nx*deltak,'FlatTime',rf_dur);
gxCrush=mr.makeTrapezoid('x',sys,'Area',-gexc.area/2);
adc = mr.makeAdc(Nx,sys,'Duration',gexc.flatTime,'Delay',gexc.riseTime);

gacqAreas=-Ny/2*deltak:deltak:Ny/2*deltak;
grewind=mr.makeTrapezoid('x',sys,'Area',(gexc.amplitude*gexc.fallTime)/2);

% - Delay
delay1 = round((TE/2-mr.calcDuration(gexc)/2-2.8e-4-mr.calcDuration(gz)/2)/seq.gradRasterTime)*seq.gradRasterTime;
delay2 = round((TE/2-mr.calcDuration(gz)/2-mr.calcDuration(gzCrush)-mr.calcDuration(gRO)/2)/seq.gradRasterTime)*seq.gradRasterTime;
delayTR = round((TR-(mr.calcDuration(gexc)+2.8e-4+mr.calcDuration(gz)+mr.calcDuration(gzCrush)+mr.calcDuration(gRO)+delay2+delay1+mr.calcDuration(gxCrush)))/seq.gradRasterTime)*seq.gradRasterTime;
% Create Blocks and loop over acqusition

for i=1:Ny
    seq.addBlock(rf,gexc)
    % seq.addBlock(gexcRe)
    seq.addBlock(mr.makeDelay(delay1))
    seq.addBlock(gzCrush)
    seq.addBlock(rfref,gz)
    phase = mr.makeTrapezoid('y',sys,'Area',gacqAreas(i),'Duration',1.5e-3);
    seq.addBlock(grewind,phase,gzCrush)
    %seq.addBlock(grewind,gzCrush)
    seq.addBlock(mr.makeDelay(delay2))
    seq.addBlock(gexc,adc)
    seq.addBlock(gxCrush)
    seq.addBlock(mr.makeDelay(delayTR))
end
% seq.addBlock(mr.makeDelay(1))

name = ['SPEN_SE_' datestr(datetime("today"))];
rep = check(seq,TE);
fprintf([rep{:}]);
% R
seq.setDefinition('Name', name);
%seq.setDefinition('FOV', fov);
seq.setDefinition('FOV', [fov, fov, sliceThickness])
seq.setDefinition('sweepBw',sweepBw);
seq.setDefinition('rf_dur',rf_dur);
seq.setDefinition('Ny',Ny);
seq.setDefinition('Nx',Nx);
seq.setDefinition('Gexc',gexc.amplitude);
seq.setDefinition('TE',TE);
seq.setDefinition('TR',TR);

seq.write([name '.seq']);
%seq.install('siemens')

R=gexc.flatArea/Nx*fov
return
%% Reconstruction for R=1
clear all
close all
clc
[scan_ID,path] = uigetfile('*.dat','select the rawdata*.dat file','C:\Users\andre\Desktop\Messungen_fuer_MA\SPEN_SE');
s = convertStringsToChars(strcat(path,scan_ID));
raw=mapVBVD(s);
rawdata=permute(raw.image{''},[3 1 2]);

reco(:,:,1:size(rawdata,3))=ifftshift(ifft2(ifftshift(rawdata(:,:,1:end))));
sos=sqrt((sum(abs(reco(:,:,1:size(reco,3))).^2,3)));

imagesc(sos);axis('tight','equal'); colormap('gray')

%% Fresnel reconstruction for R>1
clear all
close all
clc
[scan_ID,path] = uigetfile('*.dat','select the rawdata*.dat file','C:\Users\andre\Desktop\Messungen_fuer_MA\SPEN_SE');
s = convertStringsToChars(strcat(path,scan_ID));
raw=mapVBVD(s);
rawdata=permute(raw.image{''},[3 1 2]);
seq_name=[s(1:end-4) '.seq'];
sys = mr.opts('MaxGrad',40,'GradUnit','mT/m',...
    'MaxSlew',200,'SlewUnit','T/m/s',...
    'rfRingdownTime', 20e-6, 'rfDeadtime', 100e-6,'B0',3);
seq=mr.Sequence(sys);
seq.read(seq_name);
fov = seq.getDefinition('FOV'); % The FOV of phase encoded dimension
rf_dur = seq.getDefinition('rf_dur'); % The duration of the chirp pulse
Ny = seq.getDefinition('Ny'); % The number of points in the phase encoded dimension
gexc = seq.getDefinition('Gexc'); % Excitation gradient amplitude
R=gexc*rf_dur/Ny*fov; % R-factor
sweepBw=seq.getDefinition('sweepBw'); % Sweeping bandwidth

% Fresnel convolution reconstruction. Rewriting FFT as a convolution: 
% Zaitsev, M., Schultz, G., Hennig, J., Gruetter, R. and Gallichan, D. (2015),
% Parallel imaging with phase scrambling. Magn. Reson. Med., 73: 1407-1419.
% https://doi.org/10.1002/mrm.25252 
% Reconstruction algorithm can be requested: maxim.zaitsev@uniklinik-freiburg.de
for i=1:size(rawdata,3)
    reco(:,:,i)=(new_fresnel_conv(squeeze(rawdata(:,:,i)),[0,-R],[0,R],[0,0.1]));
end

image=sqrt((sum(abs(reco(:,:,1:size(reco,3))).^2,3)));
figure; imagesc(imrotate(image,180)); colormap('gray');axis('equal')
