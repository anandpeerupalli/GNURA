% [lenDoppOS,min(DoppOS),max(DoppOS), sat_el1,sate_l2] = VESTA_DopplerCalcV10( ...
%     tstart, passdurnsec, SHLatLong, GSLatLong, fwdpath,  fsband, fvhf, Doppname, demopass)
%
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

demopass = 1;

format compact
% saveFolder  = '/Users/arunasmacikunas/Documents/MATLAB/exactEarth/MT Generation';

sat                 = 37793;        % change this later to the NORAD ID for VESTA

fs                  = 9.6e3 * 4 * 8;        % will create a 500 GB file and require twice as much RAM - do not run!
fs                  = 9.6e3 * 4;
fs                  = 4;
fs                  = 40;
fs                  = 9.6e3;
% fs          = 26 * 9.6e3;           % approx. 7.4 kHz margin to min/max Doppler in Fs/4 units
% fs          = 100;

if(demopass == 1)
    vesta           = 0;        % assume EV6
    
    vesta           = 1;        % assume VESTA
    
%     fs              = 24 * 9.6e3;
    
    if(vesta)
        sat             = 43781;        % EV6 for times used below
        tt = 7.374078302008309e+05 - ( (-2*5874 + 62*60) / 86400 ); % evening pass 13 Dec 2018
        tt  =  tt + ( 60 / 86400 );      % adjust time about +70 seconds

    else
        sat             = 6;        % EV6 for times used below
        tcent =  7.371617710817399e+05 + 1.063555556e+02 + ( 20 / 86400 );         % this is the time used for good simulation for sample TLE ~6:30 pm 11 April 2018 + about 100 days for EV6
        tt  =  tcent - ( (passdurnsec/2) / 86400 );         % this is the time used for good simulation for sample TLE ~6:30 pm 11 April 2018 for EV6
        tt  =  tt + ( 10 / 86400 );         % adjust time + 10 seconds
        tt  =  tt - ( 5810 / 86400 );         % adjust time about -1 orbit, 6000 seconds
        %     tt  =  tt + ( 5810 / 86400 );   % adjust time about 1 orbit, 6000 seconds
        tt  =  tt - ( 14*5874 / 86400 );      % adjust time about -15 orbit, 6000 seconds

        
    end
    %     passdurnsec      = 12 * 60;      % 14 minutes is below horizon, 12 minutes is about 4 deg min elevation
    passdurnsec     = 9.5 * 60;      % 10 minutes duration satisfies minimum elevation of about 8 degrees, 9.5 min. about 9 deg.
    
    passdurnsec     = 8 * 60;      % 10 minutes duration satisfies minimum elevation of about 8 degrees, 9.5 min. about 9 deg.
    
    
    % Goonhilly 51.0N, 5.3W
    GSLatLong       = [51, 5.3, 0.050];             % Goonhilly Earth Station location, assume 50 m elevation
    
    % Harwich 51.947535 1.28864
    SHLatLong       = [51.9475, 1.28864, 0.050];    % Ship located near Harwich harbour, assume 50 m elevation
    Doppname        = 'DoppTest';
    
else
    tt = tstart;           % note - this is never correct, the time is the start of pass we plan to use
    fs          = 24 * 9.6e3;     % Operational standard - approx. 3 kHz margin to min/max Doppler in Fs/4 units (to DC and Fs/2)
    
end

% GSLatLong       = pi/180 * GSLatLong;
% SHLatLong       = pi/180 * SHLatLong;

demopass            = 0;

Doppname        = 'DoppTestFwd9p6Test';
% frequencies, etc. to assume for forward link
fsband      = 2234.5e6;
fvhf        = 157.3825e6;

[lenDoppOS, minDopp, maxDopp, DoppOS, sat_el1, sat_el2, mattime] = VESTA_DopplerCalcV10b( ...
    sat, tt, passdurnsec, fs, SHLatLong, GSLatLong, fsband, fvhf, Doppname, demopass);

Doppname        = 'DoppTestRtn9p6Test';
% frequencies, etc. to assume for return link

fsband      = 2060.0e6;
fvhf        = 161.8625e6;

[lenDoppOS, minDopp, maxDopp, DoppOS, sat_el1, sat_el2, mattime] = VESTA_DopplerCalcV10b( ...
    sat, tt, passdurnsec, fs, SHLatLong, GSLatLong, fsband, fvhf, Doppname, demopass);

[idx dump]  = find(sat_el1>0);
idx         = idx(1:end-1);
max(sat_el1)
max(sat_el2)

figure(100); plot(mattime,sat_el1,'.k'); grid; title('Satellite Elev (deg) (gateway K, ship R)');  ylabel('Elev (deg)'); datetick('x',13,'keepticks')
hold on; plot(mattime,sat_el2,'.r'); hold off

