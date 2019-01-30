function savewf(signal, nam, savWav, savBin, savDec, savFloat, saveFolder, backoffdB, fs, wavBitDepth, endian )
%
% copyright (c) Arunas Macikunas, Waves in Space Corp. 2018
%   Ver. 0.1    9 April 2018 Initial
%   Ver. 0.2    16 April 2018 added features for GNU output (big endian), 8
%   bit interleaved binary, etc.
%
matver      = 'R12a';

if(endian == 'b')
    floatendian = 'ieee-be';
else
    floatendian = 'ieee-le';
end

stimLenfull         = length(signal)/fs;
ver                 = 1;
saveName            = sprintf('%s_%.0fsec_fs%.0f_%sHzVer%.0f',nam,stimLenfull,fs,ver,ver);

% subfolder = '/wave files/';
wavSignal           = [real(signal),imag(signal)];

normalizeDenom      = max([max(wavSignal(:,1)),-min(wavSignal(:,1)),max(wavSignal(:,2)),-min(wavSignal(:,2))]);
dBadjustment        = - backoffdB;
scaling             = 10^(dBadjustment/20);
normalizeDenom      = normalizeDenom/scaling;

if normalizeDenom > 1
    wavSignal(:,1)  = wavSignal(:,1)/normalizeDenom;
    wavSignal(:,2)  = wavSignal(:,2)/normalizeDenom;
end

if(savWav)
    
    if ~exist([saveFolder,'/wave files/'],'dir')
        mkdir([saveFolder,'/wave files/'])
    end
    
    if(matver == 'R12a')
        wavwrite(wavSignal, fs, wavBitDepth, strcat(saveFolder,'/wave files/', saveName, '.wav'));
    else
    % audiowrite(filename, C, Fs)
        audiowrite(strcat(saveFolder,'/wave files/',saveName,'.wav'), wavSignal, fs, 'BitsPerSample', wavBitDepth)
    end
    
end

dataComplex = wavSignal(:,1) + wavSignal(:,2)*1i;
decIn = zeros(length(dataComplex)*4,1);

decIn(1:4:end-3)=real(dataComplex);
decIn(2:4:end-2)=imag(dataComplex);
decIn(3:4:end-1)=real(dataComplex);
decIn(4:4:end)=imag(dataComplex);

decIn = decIn/max(decIn);

if(savDec)
    if ~exist([saveFolder,'/DecIn/'],'dir')
        mkdir([saveFolder,'/DecIn/'])
    end
    
    decIn = decIn*2^15;
    filename = sprintf('%s/DecIn/%s.0s',saveFolder,saveName);
    fid = fopen(filename,'w');
    fwrite(fid,decIn,'int16');
    fclose(fid);
end

if(savBin||savFloat)
    
    decIn               = zeros(length(dataComplex)*2,1);
    
    decIn(1:2:end-1)    = real(dataComplex);
    decIn(2:2:end)      = imag(dataComplex);
    
    decIn               = decIn/max(decIn);
    
    if ~exist([saveFolder,'/BinaryIn/'],'dir')
        mkdir([saveFolder,'/BinaryIn/'])
    end
end

if(savBin)

    decIn8                = decIn*2^7;
    decIn16               = decIn*2^15;
    
    saveName = sprintf('%s_%.0fsec_%.0f_%s16bVer%.0f',nam,stimLenfull,fs,ver);
    filename            = sprintf('%s/BinaryIn/%s.bin',saveFolder,saveName);
    fid                 = fopen(filename,'w', endian);
    fwrite(fid,decIn16,'int16', endian);
    fclose(fid);
    
    saveName = sprintf('%s_%.0fsec_%.0f_%s8bVer%.0f',nam,stimLenfull,fs,ver);
    filename            = sprintf('%s/BinaryIn/%s.bin',saveFolder,saveName);
    fid                 = fopen(filename,'w', endian);
    fwrite(fid,decIn8,'int8', endian);
    fclose(fid);
end

if(savFloat)
    if ~exist([saveFolder,'/float files/'],'dir')
        mkdir([saveFolder,'/float files/'])
    end
    
    saveName        = sprintf('%s_%.0fsec_fs%.0f_%sHzVer%.0f32',nam,stimLenfull,fs,ver,ver);
    filename        = sprintf('%s/float files/%s.flt',saveFolder,saveName);
    fid             = fopen(filename,'w');
%     fwrite(fid,wavSignal,'float32',floatendian);
%     fwrite(fid,wavSignal,'float',floatendian);
    fwrite(fid,decIn*scaling,'float',floatendian);

    fclose(fid);
    
    % there is no complex 64 file format, only variable format
%     saveName        = sprintf('%s_%.0fsec_fs%.0f_%sHzVer%.0f64',nam,stimLenfull,fs,ver,ver);
%     filename        = sprintf('%s/float files/%s.flt',saveFolder,saveName);
%     fid             = fopen(filename,'w');
%     fwrite(fid, signal/normalizeDenom, 'complex64');
%     fclose(fid);
end
