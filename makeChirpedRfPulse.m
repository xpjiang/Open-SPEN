%% Open SPEN in Pulseq - Chirp RF pulse function 
% Andreas Holl
% Division of Medical Physics, Department of Diagnostic and Interventional Radiology,
% University Medical Center Freiburg, Faculty of Medicine, University of Freiburg, Freiburg, Germany
% Email: andreas.holl@uniklinik-freiburg.de
% March. 23, 2024

% A part of the following code was taken from 
%https://github.com/mikgroup/sigpy/tree/main/sigpy/mri.

function [rf, pm, gz, gzr, delay] = makeChirpedRfPulse(varargin)

validPulseUses = mr.getSupportedRfUse();
persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'makeChirpedRfPulse';
   
    % RF params
    % addRequired(parser, 'type', @(x) any(validatestring(x,validPulseTypes)));
    addOptional(parser, 'system', mr.opts(), @isstruct);
    addParamValue(parser, 'duration', 0, @isnumeric);
    addParamValue(parser, 'ang', 0, @isnumeric);
    addParamValue(parser, 'freqOffset', 0, @isnumeric);
    addParamValue(parser, 'phaseOffset', 0, @isnumeric);
    addParamValue(parser, 'beta', 800, @isnumeric);
    addParamValue(parser, 'mu', 4.9, @isnumeric);
    addParamValue(parser, 'n_fac', 40, @isnumeric);
    addParamValue(parser, 'bandwidth', 0, @isnumeric);
    addParamValue(parser, 'adiabaticity', 0.4, @isnumeric);
    % Slice params
    addParamValue(parser, 'maxGrad', 0, @isnumeric);
    addParamValue(parser, 'maxSlew', 0, @isnumeric);
    addParamValue(parser, 'sliceThickness', 0, @isnumeric);
    addParamValue(parser, 'delay', 0, @isnumeric);
    addParamValue(parser, 'dwell', 0, @isnumeric); % dummy default value
    % whether it is a refocusing pulse (for k-space calculation)
    addOptional(parser, 'use', '', @(x) any(validatestring(x,validPulseUses)));
    % addParamValue(parser, 'freqPPM', 0, @isnumeric);
    % addParamValue(parser, 'phasePPM', 0, @isnumeric);
end
parse(parser, varargin{:});
opt = parser.Results;
opt.duration=opt.duration;%-1e-5
if opt.dwell==0
    opt.dwell=opt.system.rfRasterTime;
end
Nraw = round(opt.duration/opt.dwell+eps);
N = floor(Nraw/4)*4; % number of points must be divisible by four -- this is a requirement of the underlying library
t=(0:N-1)*opt.duration/N;
am=1-abs(cos(pi*t/opt.duration)).^opt.n_fac;
% am=ones(1,length(am));
fm=linspace(-opt.bandwidth/2,opt.bandwidth/2,N)*2*pi;
pm=cumsum(fm)*opt.dwell;
ifm=length(fm)/2;
dfm=fm(length(fm)/2);
%[dfm,ifm]=min(abs(fm)); % find the center of the pulse
% we will also use the ocasion to find the rate of change of the frequency
% at the center of the pulse
if dfm==0
    pm0=pm(ifm);
    am0=am(ifm);
    roc_fm0=abs(fm(ifm+1)-fm(ifm-1))/2/opt.dwell;
else
    % we need to bracket the zero-crossing
    if fm(ifm)*fm(ifm+1) < 0
        b=1;
    else
        b=-1;
    end
    pm0=(pm(ifm)*fm(ifm+b)-pm(ifm+b)*fm(ifm))/(fm(ifm+b)-fm(ifm));
    am0=(am(ifm)*fm(ifm+b)-am(ifm+b)*fm(ifm))/(fm(ifm+b)-fm(ifm));
    roc_fm0=abs(fm(ifm)-fm(ifm+b))/opt.dwell;
end
pm=pm-pm0;
a=((roc_fm0*opt.adiabaticity)^0.5/2/pi/am0);
signal = ((a*am.*exp(1i*pm))*1.051/90)*opt.ang;
if (N~=Nraw)
    % we need to pad the signal vector
    Npad=Nraw-N;
    signal=[zeros(1,Npad-floor(Npad/2)) signal zeros(1,floor(Npad/2))];
    N=Nraw;
end
%BW = opt.timeBwProduct/opt.duration;
t = ((1:N)-0.5)*opt.dwell;
%flip = abs(sum(signal))*opt.dwell*2*pi;
rf.type = 'rf';
rf.signal = signal;
rf.t = t;
rf.shape_dur=N*opt.dwell;
rf.freqOffset = opt.freqOffset;
rf.phaseOffset = opt.phaseOffset;
rf.deadTime = opt.system.rfDeadTime;
rf.ringdownTime = opt.system.rfRingdownTime;
rf.delay = opt.delay;


% rf.center = rf.shape_dur / 2;
% 
% rf.freqPPM = opt.freqPPM;
% rf.phasePPM = opt.phasePPM;


if ~isempty(opt.use)
    rf.use=opt.use;
else
    rf.use='excitation';
end
if rf.deadTime > rf.delay
    rf.delay = rf.deadTime;
end
if rf.ringdownTime > 0 && nargout > 1
    delay=mr.makeDelay(mr.calcDuration(rf)+rf.ringdownTime);
end

