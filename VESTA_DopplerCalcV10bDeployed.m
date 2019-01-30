function [lenDoppOS, minDopp, maxDopp, DoppOS, sat_el1, sat_el2, mattime] = VESTA_DopplerCalcV10b( ...
    sat, tstart, passdurnsec, fs, SHLatLong, GSLatLong,  fsband, fvhf, Doppname, demopass)

% [lenDoppOS, minDopp, maxDopp, DoppOS, sat_el1, sat_el2, mattime] = VESTA_DopplerCalcV10b( ...
%     sat, tt, passdurnsec, fs, SHLatLong, GSLatLong, fsband, fvhf, Doppname, demopass);%
% Frequency correction waveform generation for VESTA
% Arunas Macikunas, Copyright (c) 2018 Waves in Space Corporation
% Ver. 0.0 11 April 2018
% Ver. 0.1 added dual station total Doppler
% Ver. 0.2 added CW waveform generation, and efficiencies for large arrays
% (28 April 2018)
% Ver. 0.3 added formulas for forward and return total Doppler for ship and
%          gateway relevant frequencies
% Ver. 0.4 changed sampling rates and variables for Doppler vs. CW time file, and CW
%          correction output vectors
% Ver. 0.7 main signal generation at offset frequency (fs/2 offset) and using VCO function
% Ver. 0.8 (V8) change to start time and duration notation for call statement and frequency calculation vector method using
%      integrated phase rather than frequency, *** Dopp freq scaled back by 40X, and Fs is too low - illustration only ***
% Ver. 0.9 (V9) final version, fs = 26 X 9.6 kS/s, full Doppler, shorter pass, signal1 deleted, final script version
% Ver. 1.0 (V10) function version
%
% Function statement example syntax below
%    in functional form use [outputarguments] = VESTA_DOPP(tstart, passdurnsec, SHLatLong, GSLatLong, fwdpath, wavfilenam, demopass, ...)
%    Units lat, long in degrees and decimals, frequencies in Hz
%
% sgp4init

% demopass = 1;

format compact
saveFolder  = '/Users/arunasmacikunas/Documents/MATLAB/exactEarth/MT Generation';

debugflg        = 1;        % = 1, debug mode - do not write output files
debugflg        = 0;

savWav          = 1;
savBin          = 0;
savDec          = 0;
savFloat        = 0;
% savConj         = 0;

if(demopass == 1)
    % frequencies, etc. to assume
    fwdpath         = 1;        % select to calculate
    if(fwdpath==1)
        fsband      = 2234.5e6;
        fvhf        = 157.3825e6;
    else
        fsband      = 2060.0e6;
        fvhf        = 161.8625e6;
    end
    
    %     passdurnsec      = 12 * 60;      % 14 minutes is below horizon, 12 minutes is about 4 deg min elevation
    passdurnsec      = 9.5 * 60;      % 10 minutes duration satisfies minimum elevation of about 8 degrees, 9.5 min. about 9 deg.
    
    passdurnsec      = 8 * 60;      % 10 minutes duration satisfies minimum elevation of about 8 degrees, 9.5 min. about 9 deg.
    
    tcent =  7.371617710817399e+05 + 1.063555556e+02 + ( 20 / 86400 );         % this is the time used for good simulation for sample TLE ~6:30 pm 11 April 2018 + about 100 days for EV6
    tt  =  tcent - ( (passdurnsec/2) / 86400 );         % this is the time used for good simulation for sample TLE ~6:30 pm 11 April 2018 for EV6
    tt  =  tt + ( 10 / 86400 );         % adjust time + 10 seconds
    tt  =  tt - ( 5810 / 86400 );         % adjust time about -1 orbit, 6000 seconds
    %     tt  =  tt + ( 5810 / 86400 );         % adjust time about 1 orbit, 6000 seconds
    tt  =  tt - ( 14*5874 / 86400 );         % adjust time about -15 orbit, 6000 seconds
    
    % Goonhilly 51.0N, 5.3W
    GSLatLong       = [51, 5.3, 0.050];             % Goonhilly Earth Station location, assume 50 m elevation
    
    % Harwich 51.947535 1.28864
    SHLatLong       = [51.9475, 1.28864, 0.050];    % Ship located near Harwich harbour, assume 50 m elevation
    Doppname        = 'DoppTest';
    
else
    tt = tstart;           % note - this is never correct, the time is the start of pass we plan to use
end

GSLatLong       = pi/180 * GSLatLong;
SHLatLong       = pi/180 * SHLatLong;

%tt = tt + 0.1/(24*3600);    % add 0.1 second - not sure why - this was bypassed Dec. 4, 2018

% GSLatLon      = [pi/180*54,pi/180*6,500];
% GSLatLon      = [pi/180*54,pi/180*6,-100];

% [sat_el,sat_az,pos,range,rangeRate] = satElevation(tt+(1:256*1440)/(1440/.5),6,GSLatLon);
% [sat_el,sat_az,pos,range,rangeRate] = satElevation(tt+(306300:306310)/(1440/.5),6,GSLatLon);
% [sat_el,sat_az,pos,range,rangeRate] = satElevation(tt+(306303:306305)/(2*1440),6,GSLatLon);

%
% set the sampling rate for calls to SGP4 Propagation, and Doppler
% prediction function created
%
% sampspersec     = 100;        % quite a lot of 'jitter' in the results
% sampspersec     = 10;         % lot of 'jitter' in the results
% sampspersec     = 0.25;       % minimal visisble 'jitter' in the results
sampspersec     = 1;            % not much jitter 'jitter' in the results

% timevec         = (-passdurnsec/2*sampspersec:passdurnsec/2*sampspersec)/(1440*60*sampspersec);

timevec         = (0:passdurnsec*sampspersec) / (1440*60*sampspersec);

mattime         = tt + timevec;
lenmattime      = length(mattime);

% mattimeextra    = mattime(end) + 1 /(1440*60*sampspersec);

[a b] = find( abs(diff(diff(mattime))>0));

% Note:  ***** zero out points in mattime that have 'bumps' - but this will require
% changes to interpolation based on time jumps
% 2nd idea - try to replace these values with interpolated time from
% adjacent time samples *****

% occassional errors due to round-off are of the order of +/- 1.15e-10
% every 17 or 18 samples, this error is equivalent to 1.15e-10 * 86400
% (seconds/day) or about 10 microseconds - jumps in TLE appear to be
% equivalent to about 1 second error in time - so the errors in range rate
% must be either internal to the SGP4 propagator (likely due to quantization
% effects etc. amplifying these small inconsistencies in time base)

if( 1==0)
    
    %     load VESTA_Doppler_variables_4Dec2018
    load VESTA_Doppler_variables_9Dec2018
    
else
    
    [sat_el1,sat_az1,pos1,range1,rangeRate1] = satElevation(mattime,sat,GSLatLong);     % gateway station S-band
    [sat_el2,sat_az2,pos2,range2,rangeRate2] = satElevation(mattime,sat,SHLatLong);     % ship station VHF band
    
end

% save VESTA_Doppler_variables_9Dec2018

% note:
% rangerate, range, etc. 1 = Gateway to/from satellite (i.e. S-band link)
% rangerate, range, etc. 2 = Ship to/from satellite (i.e. VHF link)
%

% os_factor   = 1;
% fitcoefph1		    = polyfit(1:lenmattime, range1', 9);		    	% linear array only, fit to coef. = 5rd order
% range1F			    = polyval(fitcoefph1,(os_factor:lenmattime*os_factor)/os_factor);
range1F          = medfilt1(range1,5);      % median filter range function to remove some bumps

if(1==0)
    figure(95); plot(mattime(1:(end-1)),diff(range1)./diff(timevec),'.k'); grid; title('Satellite Elev Range Rate and Calculated)');  ylabel('Range Rate (km/s)'); datetick('x',13,'keepticks')
    hold on; plot([mattime(1:(end-1))],rangeRate1,'.r');        % testing shows that range rate is smooth, delta(R)/delta(T) has inconsistency in the middle
end

% plot(mattime(1:(end-1)),diff(range1F)/(mean(diff(timevec))),'.m');
plot(mattime(1:(end-1)),diff(range1F)/(mean(diff(timevec))),'.g'); hold off

rangeRate1       = rangeRate1/(1440*60);  % change units to km/s
rangeRate2       = rangeRate2/(1440*60);  % change units to km/s

if(1==0)
    rangeRate1      = medfilt1(rangeRate1,9);
    rangeRate2      = medfilt1(rangeRate2,9);
end
if(1==0)
    rangeRate1      = changefilter(rangeRate1,1.75,0.5);
    rangeRate2      = changefilter(rangeRate2,1.75,0.5);
end

% rangeRate1      = medfilt1(rangeRate1,11);
% rangeRate2      = medfilt1(rangeRate2,11);
%
% [sat_el,sat_az,pos,range,rangeRate] = satElevation(tt+76576/720-100+(1:1440)/(1440/2),6,GSLatLon);

[idx dump]  = find(sat_el1>0);
idx         = idx(1:end-1);
max(sat_el1)
max(sat_el2)

figure(100); plot(mattime,sat_el1,'.k'); grid; title('Satellite Elev (deg) (gateway K, ship R)');  ylabel('Elev (deg)'); datetick('x',13,'keepticks')
hold on; plot(mattime,sat_el2,'.r'); hold off

figure(101); plot(range1,'.k'); grid; title('Range (km, elev > 0 is red)')
hold on; plot(idx, range1(idx),'.r'); plot(idx, range2(idx),'.g'); hold off

figure(102); plot(rangeRate1,'.k'); grid; title('Range rate (km/s, elev > 0 is red)');ylabel('km/s');
hold on; plot(idx, rangeRate1(idx),'.r');  plot(rangeRate2,'.g');  hold off

DoppV     = 1e3*rangeRate2/3e8*fvhf;       % VHF Doppler in Hz (1e3 changes range rate scaling in km/s to m/s)
doppnorm    = DoppV/max(abs(DoppV));
DoppS     = 1e3*rangeRate1/3e8*fsband;

if(1==0)
    figure(103); plot(DoppV/1e3,'.k'); grid; title('VHF Doppler (kHz, elev > 0 is red)');ylabel('Frequency Offset (kHz)');
    hold on; plot(idx, DoppV(idx)/1e3 ,'.r'); hold off
    
    figure(104); plot(DoppS/1e3,'.k'); grid; title('S-band Doppler (kHz, elev > 0 is red)');ylabel('Frequency Offset (kHz)');
    hold on; plot(idx, DoppS(idx)/1e3 ,'.r'); hold off
end

figure(105); plot( (DoppV+DoppS)/1e3,'.k'); grid; title('Total Doppler (kHz, elev > 0 is red)');ylabel('Frequency Offset (kHz)');
hold on; plot(idx, (DoppV(idx)+DoppS(idx))/1e3 ,'.r'); hold off

% fs          = sampspersec;
backoffdB   = 1;
wavBitDepth = 16;
endian      = 'l';      % endian - either l or b, for big or little endian for 2 byte binary words, Mac is little, Unix is big

DoppT       = (DoppV + DoppS);
% DoppT       = (0.5 * DoppV + 0*DoppS);

lenDopp     = length(DoppT);   % 1 second long

if(1==1)
    %     Doppfilt = changefilter(DoppT,3.5,0.5);
    %     Doppfilt = changefilter(DoppT,1.5,0.5);
    %     Doppfilt = changefilter2(DoppT,1.5,1);
    %     Doppfilt = changefilter3(DoppT,.5,1.25);
    
    Doppfilt = changefilter3(DoppT,1.5,2);
end

DoppT       = Doppfilt;

if(1==0)
    os_factor   = 1;
    fitcoefph1		    =	polyfit(1:lenDopp, DoppT'-min(DoppT), 7);		    	% linear array only, fit to coef. = 5rd order
    DoppTS			    =	min(DoppT) + polyval(fitcoefph1,(os_factor:lenDopp*os_factor)/os_factor);
    
    % figure(106); plot(DoppT,'.'); grid; title('DoppT and smoothed');hold on; plot(DoppTS,'.r'); plot(10*(DoppT-DoppTS),'.g'); hold off
    figure(106); plot(DoppT,'.'); grid; title('DoppT and smoothed');hold on; plot(DoppTS,'.r'); hold off
end

signal1     = DoppT / max(abs(DoppT));

signal1r    = -signal1;

if(sampspersec >= 0.01 && ~debugflg)
    nam         = 'Doppler3_total_inkHz';
    
    savewfR12a(DoppT/1e3, nam, savWav, savBin, savDec, savFloat, saveFolder, backoffdB, sampspersec, wavBitDepth, endian)
    nam         = 'Doppler3_total_norm';
    
    savewfR12a(signal1, nam, savWav, savBin, savDec, savFloat, saveFolder, backoffdB, sampspersec, wavBitDepth, endian)
    nam         = 'Doppler3_total_Rnorm';
    
    savewfR12a(signal1r, nam, savWav, savBin, savDec, savFloat, saveFolder, backoffdB, sampspersec, wavBitDepth, endian)
end

% create high-sample rate Doppler correction file (complex sine wave)

% create real-time SIN wave correction function
% periodsin2  = fs/CWfreq2;
% signal      = exp(-(2*1j*pi*((1:lensig)-1))/periodsin2).';

interprate  = fs / sampspersec;         % rate of oversampling of CW correction function over Doppler function
CWfreq      = 1e3;

% lensig      = fs;   % 1 second long
% freqsin     = CWfreq / fs;
% signal      = exp(-(2*1j*pi*((1:lensig)-1)).*freqsin).';

freqsin     = CWfreq / fs;

fprintf('Basic sampling rate of orbit propagator %.2f (Hz), pass length %.1f (sec), desired rate %.1f (kHz), overampling rate %.0f \n', ...
    sampspersec, lenDopp/sampspersec, fs/1e3, interprate);

if(1==0)        % VARIABLES BELOW NOT USED
    time1       = 1:length(DoppT);
    timeinter   = (1:length(DoppT)*interprate)/interprate;
end

% Note:  all approaches below have serious artifacts (100s of Hz + of fast
%       ripple type effects in interplotated files0
% ts1 = resample(ts, time, interp_method) ;
% p           = pchip(x,y,xq)
% signal1     = resample(DoppT, 1, 4) ;
% signal1      = pchip(time1,DoppT,timeinter);
%
lenOS       = (lenDopp-1)*interprate;
DoppOS      = ones(1,lenOS);

DoppOS(lenOS) = -1;

% ************** Reduce Doppler 20 fold for Testing only at low sampling rates ***************
% DoppT = DoppT / 20;
% ************** Reduce Doppler 20 fold for Testing only at low sampling rates ***************


for ii = 1:lenDopp-1
    
    subidx  = 1:(interprate);
    idx     = (ii - 1) * interprate + subidx;
    
    DoppSlope = ( DoppT(ii+1) - DoppT(ii) )/interprate;
    %     DoppSlope = ( DoppT(ii+1) - DoppT(ii) );
    
    DoppOS(idx) = DoppT(ii) + (subidx - 1) * DoppSlope;
    
    %     if ii == 112
    %        figure(116); plot((diff(DoppOS(1:min(1e7,lenOS)))),'.');grid; title('Diff DoppOS (= Doppler offset) Hz');
    %
    %        fprintf('look here!/n');
    %     end
    %     for ij = 1:interprate
    %         indx = indx + 1;
    %         DoppOS(indx) = DoppT(ii) + ( DoppT(ii+1) - DoppT(ii) )/interprate;
    %     end
end

% sine was is
DoppOS(lenOS)   = DoppT(end);
lenDoppOS       = length( DoppOS );

% DoppOS          = DoppOS / fs;
DoppOS          = DoppOS;

% signal1      = exp(-(2*1j*pi*((1:length(DoppOS))-1)) .* DoppOS/fs).';

% signal1      = exp( -(2*1j*pi*(DoppOS .* (1:length(DoppOS))-1)) / ( length(DoppOS) *fs ) ).';   % updated Dec. 4, 2018

% DoppOS          = 5e2 + 200 * sin( 2*pi*500*((1:(lenDoppOS-1))/lenDoppOS));
% DoppOS          = 44 + 11 * sin( 2*pi*10*( (1:(lenDoppOS-1)) /lenDoppOS ) );

% Notes:

% 4 Dec. 2018
% something wrong with argument - check scaling for CW - do FFT, etc.,
% shorter sequence, known sampling rate
% continuously changing frequency should be fine, try FSK - two frequencies
% also, see if Matlab has a 'VCO' function, this should produce CW from
% control voltage that is proportional to desired frequency
%

% DoppOS          = 1000 + 1000 * sin( 2*pi*10*( (1:(lenDoppOS)) /lenDoppOS ) );

minDopp         = min(DoppOS) , maxDopp         = max(DoppOS)

DoppOS(1)/1e3, DoppOS(end)/1e3

if(1==1)
    
    % signal1         = exp( 1j*(2*pi*((DoppOS + min(DoppOS)) .* (1:(lenDoppOS)))) / (  2*fs ) ).';   % updated Dec. 4, 2018 - still not right, not 2-sided
    if(1==1)
        signal2         = vco(DoppOS/max(maxDopp,-minDopp), [minDopp maxDopp] + fs/4, fs);      % fs/4 offset CW wave
    else
        signal2         = vco(DoppOS/max(maxDopp,-minDopp), [minDopp maxDopp] + 0*fs/4, fs);    % good for audio test by ear
    end
    
end

% signal2b            = signal2;
% save temp1 signal2b
%

figure(91);
spectrogram(signal2,hamming(256),0,1024,fs,'yaxis')

% signal1         = vco( sawtooth(2*pi*

% if(lenOS >= 100e6)        % don't plot if there are many samples
%     clear DoppOS
% end

% nam         = 'VHFpStotal_OS';
savFloat    = 0;

if(~debugflg)
    savewfR12a(signal2, Doppname, savWav, savBin, savDec, savFloat, saveFolder, backoffdB, fs, wavBitDepth, endian)
end

if(lenOS < 10e6)        % don't plot if there are many samples
    %     figure(110); [pxx,f] = periodogram(signal1); plot(f, pxx, '.'); title('Periodogram of signal');
    %     figure(111); spectrogram(signal1,blackman(256),128,1024,fs,'yaxis'); title('Spectrogram plot of signal 1');
    figure(112); spectrogram(signal2,blackman(256),128,1024,fs,'yaxis'); title('Spectrogram plot of signal 2');
    
    figure(113); plot(DoppOS,'.');  title('DoppOS'); grid;
    figure(114); plot(real(signal1),'.');  title('Real signal1');
    
    clear pxx f
else
    clear signal1
end

if(lenOS < 10e6 && exist('DoppOS','var'))
    figure(115); plot(diff(DoppOS(1:min(1e7,lenOS))),'.');grid; title('Diff DoppOS (= Doppler offset) Hz');
    figure(116); plot(diff(diff(DoppOS(1:min(1e7,lenOS)))),'.');grid; title('Diff Diff DoppOS (= Doppler offset) Hz');
end
