%% Open SPEN in Pulseq - SPEN-SE-EPI
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
fov = 300e-3;%256e-3
Nx = 100;
Ny = 100;
deltak = 1/fov; % Pulseq toolbox defaults to k-space units of m^-1
kWidth = Nx*deltak;

% Define sequence parameters
TE = 68e-3;

% set system limits
sys = mr.opts('MaxGrad',40,'GradUnit','mT/m',...
    'MaxSlew',200,'SlewUnit','T/m/s',...
    'rfRingdownTime', 20e-6, 'rfDeadtime', 100e-6,'B0',3);

% Create a new sequence object
seq = mr.Sequence(sys);

% Calculate SPEN-conditions: 
% sweepBw*rf_dur = gexc.amplitude*rf_dur*fov = gacq.amplitude*Ny*mr.calcDuration(gacq)*fov

rf_dur = 4e-3;
sweepBw=50000;

[rf,pm] = makeChirpedRfPulse('duration',rf_dur,'delay',sys.rfDeadTime+3e-5,'bandwidth',sweepBw, ...
    'ang',90,'n_fac',40,'system',sys);
[bw,f0,M_xy_sta,F1]=mr.calcRfBandwidth(rf);

[rf,pm] = makeChirpedRfPulse('duration',rf_dur,'delay',sys.rfDeadTime+3e-5,'bandwidth',sweepBw*sweepBw/bw, ...
    'ang',90,'n_fac',40,'system',sys);
[bw,f0,M_xy_sta,F1]=mr.calcRfBandwidth(rf);

% Create SPEN-gradient
gexc = mr.makeTrapezoid('y',sys,'Amplitude',bw/fov,'FlatTime',rf_dur,'Delay',sys.rfDeadTime);

% Create a slice selective sinc 180° RF pulse and crusher
dur_ref = 3e-3;
sliceThickness = 5e-3;
[rfref, gz] = mr.makeSincPulse(pi,sys,'Duration',dur_ref,...
    'SliceThickness',sliceThickness,'apodization',0.2,'timeBwProduct',4,'PhaseOffset',0,'use','refocusing');
gzCrush = mr.makeTrapezoid('z', 'Amplitude', -gz.amplitude*3, 'Duration', (dur_ref)/2);

% Create other gradients and events
gRO = mr.makeTrapezoid('x',sys,'FlatArea',kWidth,'FlatTime',3.2e-4);%3.2e-4
adc = mr.makeAdc(Nx,sys,'Duration',gRO.flatTime,'Delay',gRO.riseTime);

    % - Phase blip
gacq = mr.makeTrapezoid('y',sys,'Area',gexc.flatArea/Ny);

    % - Rephase
gRe = mr.makeTrapezoid('y',sys,'Area',(gexc.area-gexc.flatArea)/2);

    % - Crusher
gxCrush1 = mr.makeTrapezoid('x',sys, 'Area', gRO.area/4);
gxCrush2 = mr.makeTrapezoid('x',sys, 'Area', -gRO.area/4);

    % - Delay
delay1 = round((TE/2-mr.calcDuration(gexc)/2-mr.calcDuration(gzCrush)-mr.calcDuration(rfref)/2)/seq.gradRasterTime)*seq.gradRasterTime;
delay2 = round((TE/2-mr.calcDuration(rfref)/2-mr.calcDuration(gxCrush2)-(Ny/2)*(mr.calcDuration(gRO)+mr.calcDuration(gacq))+0.5*mr.calcDuration(gRO))/seq.gradRasterTime)*seq.gradRasterTime;

% Create Blocks and loop over acqusition
seq.addBlock(rf,gexc)
seq.addBlock(mr.makeDelay(delay1))
seq.addBlock(gzCrush,gxCrush1)
seq.addBlock(rfref,gz)
seq.addBlock(gzCrush,gxCrush2,gRe)
seq.addBlock(mr.makeDelay(delay2))

for i=1:Ny
    seq.addBlock(gRO,adc);            % Read one line of k-space
    seq.addBlock(gacq);
    gRO.amplitude = -gRO.amplitude;   % Reverse polarity of read gradient
end
% seq.addBlock(mr.makeDelay(1))
%% check

%Check whether the SPEN conditions are fulfilled
assert((rf_dur*sweepBw-gexc.amplitude*rf_dur*fov<1e-10)&&(gexc.flatArea*fov-gacq.area*Ny*fov<1e-10), ...
    'SPEN conditions are not fulfilled:\nrf_dur*sweepBw = %d\ngexc.amplitude*rf_dur*fov = %d\ngacq.amplitude*Ny*d*fov = %d',sweepBw*rf_dur, ...
    gexc.amplitude*rf_dur*fov,gacq.area*Ny*fov)
R=gacq.area*fov;%Reduction factor. Robustness to field inhom. and chemical shifts decreases with R
acqDur=(mr.calcDuration(gRO)+mr.calcDuration(gacq))*Ny

name = ['SPEN_EPI_SE_' datestr(datetime("today"))];
rep = check(seq,TE);
fprintf([rep{:}]);
R
seq.setDefinition('Name', name);
seq.setDefinition('FOV', fov);
seq.setDefinition('sweepBw',sweepBw);
seq.setDefinition('rf_dur',rf_dur);
seq.setDefinition('Ny',Ny);
seq.setDefinition('Nx',Nx);
seq.setDefinition('Gexc',gexc.amplitude);
seq.setDefinition('TE',TE);
seq.setDefinition('acqDur',acqDur);
seq.setDefinition('R',R);
seq.setDefinition('SeqDur',seq.duration());

seq.write([name '.seq']);
seq.install('siemens')

return

%% Reco SPEN_SE_EPI Fresnel for low R-factors
clear all
close all
clc
[scan_ID,path] = uigetfile('*.dat','select the rawdata*.dat file','C:\Users\andre\Desktop\brainSPEN\SPEN');
s = convertStringsToChars(strcat(path,scan_ID));
clearvars -except s 
raw=mapVBVD(s); % mapVBVD for reading Siemens raw-data
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
R=seq.getDefinition('R'); % R-factor
sweepBw=seq.getDefinition('sweepBw'); % Sweeping bandwidth
rawdata=permute(raw.image{''},[3 1 2]);

% Flip even k-space lines and manual gradient imperfection correction
for k=1:size(rawdata,3)
for i=2:2:length(rawdata)
    rawdata(i,:,k)=fliplr(rawdata(i,:,k));
    % shift=1;
    % rawdata(i,shift+1:end,k)=rawdata(i,1:end-shift,k);
    rawdata(i,:,k)=circshift(rawdata(i,:,k),1);
    rawdata(i-1,:,k)=circshift(rawdata(i-1,:,k),-1);

end
end

% Fresnel convolution reconstruction. Rewriting FFT as a convolution: 
% Zaitsev, M., Schultz, G., Hennig, J., Gruetter, R. and Gallichan, D. (2015),
% Parallel imaging with phase scrambling. Magn. Reson. Med., 73: 1407-1419.
% https://doi.org/10.1002/mrm.25252 
% Reconstruction algorithm can be requested: maxim.zaitsev@uniklinik-freiburg.de
for i=1:size(rawdata,3)
    reco(:,:,i)=(new_fresnel_conv(squeeze(rawdata(:,:,i)),[-R,0],1.5*[R,0],[0.1,0]));
end
image=sqrt((sum(abs(reco(:,:,1:size(reco,3))).^2,3)));
% figure; imagesc(image,[0 1.3*10^-6]); colormap('gray');axis('equal','tight')
% title('SPEN-SE-EPI',Interpreter='latex')

% Optional step. Increasing resolution of SPEN image
test=padarray(fftshift(fft2(fftshift(image))),[50 50]);
test2=ifftshift(ifft2(ifftshift(test)));
figure; imagesc(abs(test2)); colormap('gray');axis('equal','tight')

% Signal windowing
%figure; imagesc(abs(test2),[0.2 4.5]*10^-7); colormap('gray');axis('equal','tight')

title('SPEN-SE-EPI',Interpreter='latex')

%% REconstruction for very high R-factors (R>60)
clear all
close all
clc
[scan_ID,path] = uigetfile('*.dat','select the rawdata*.dat file','C:\Users\andre\Desktop\MA_MessungenNeu');
s = convertStringsToChars(strcat(path,scan_ID));

raw=mapVBVD(s); % mapVBVD for reading in Siemens raw-data
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
Nx=Ny;
gexc = seq.getDefinition('Gexc'); % Ecitation gradient amplitude
R=seq.getDefinition('R'); % R-factor
sweepBw=seq.getDefinition('sweepBw'); % RF sweeping bandwidth
totalDur=seq.getDefinition('TotalDuration'); % Total sequence duration

% Flip even k-space lines and manual gradient imperfection correction
for k=1:size(rawdata,3)
for i=2:2:length(rawdata)
    rawdata(i,:,k)=fliplr(rawdata(i,:,k));
    % shift=1;
    % rawdata(i,shift+1:end,k)=rawdata(i,1:end-shift,k);
    rawdata(i,:,k)=circshift(rawdata(i,:,k),1);
    rawdata(i-1,:,k)=circshift(rawdata(i-1,:,k),-1);

end
end

% Exponential weighting to correct for signal attenuation due to relaxation
% processes and the need for long RF pulse durations
t=linspace(0,totalDur,Ny*Nx);
relax(1:length(t))=1./exp(-t(1:end)/0.1);
for i=1:size(rawdata,3)
newraw(i,:)=reshape(rawdata(:,:,i).',1,[]);
newraw(i,:)=newraw(i,:).*relax;
newraw2(:,:,i)=reshape(newraw(i,:), size(rawdata, 2), []).';
end

blurr(:,:,1:size(rawdata,3))=ifftshift(ifft(ifftshift(newraw2(:,:,1:end),2),[],2),2); % Blurred SPEN-image

sos=sqrt((sum(abs(blurr(:,:,1:size(blurr,3))).^2,3))); % Sum of squares from SPEN coil images

% Increasing of the resolution
padsos=padarray(fftshift(fft2(fftshift(sos))),[60 50]); 

sr=ifftshift(ifft2(ifftshift(padsos)));
srcropped=sr(5:205,:);

figure;imagesc(fliplr(abs(srcropped)));axis('tight','equal');colormap('gray')
