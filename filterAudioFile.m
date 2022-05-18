%  xf = FILTERAUDIOFILE(filePath,channel,b,a,varargin)
%
%  DESCRIPTION
%  Reads a segment of the specified CHANNEL in the audio file located in 
%  FILEPATH and filters it using the filter with coefficients (B,A). In 
%  FIR filters, A = 1. The segment starts at sample SAMPLES(1) and contains
%  at total number of samples equal to SAMPLES(2). The function supports 
%  audio files in WAV and RAW formats.
%
%  The selected audio segment is split into windows of manageable size.
%  FILTERAUDIOFILE uses the overlap-save method to combine the filtered
%  responses from consecutive windows. For FIR filters (B,1), the function
%  computes the Fast-Fourier Transform (FFT) of the filter impulse response 
%  B and windowed signals (chunks), multiplies both and calculates the 
%  Inverse FFT (IFFT) to obtain the filtered waveform of consecutive chunks, 
%  combined  by overlap-save to produce the filtered segment XF. For IIR 
%  filter (B,A), the function filters the windowed signals directly in the
%  time domain using the FILTER function. The FFT method (FIR) is more
%  computationally efficient. Note also that the FILTER method (IIR) may
%  produce artifacts if (B,A) represent a very selective filter (high
%  order, low cutoff frequency).
%
%  INPUT VARIABLES
%  - filePath [string]: absolute path of audio file
%  - channel [integer number]: index of the audio channel to be filtered. 
%    Use CHANNEL = [] to filter all channels.
%  - samples [firstSample nSamples]: section of the audio file to be
%    read. The first element is the starting audio sample (firstSample =
%    1,2,...) and the second element the number of samples to be read.
%  - b [numeric vector]: filter coefficients (numerator).
%  - a [numeric vector]: filter coefficients (denominator). A = 1 for
%    FIR filters (i.e. moving-average filter) and a vector for IIR filters.
%  - sampleRate (varargin{1}) [integer number]: sampling rate (RAW).
%  - numChannels (varargin{2}) [integer number]: number of channels (RAW).
%  - bitsPerSample (varargin{3}) [integer number]: bit resolution (RAW).
%  - byteOrder (varargin{4}) [integer number]: endianess (RAW).
%
%  OUTPUT VARIABLES
%  - xf [numeric vector]: vector of filtered audio samples. The data is
%    returned as 'single' class.
%
%  INTERNALLY CALLED FUNCTIONS
%  - readwavHeader
%  - readAudioFile
%
%  CONSIDERATIONS & LIMITATIONS
%  - Use filterAudioFile instead of #filter for file sizes of the order
%    of tens of MegaBytes. The speed will be improved and the amount of
%    used memory drastically reduced.
%  - Only wav and raw formats are supported.
%  - The processing speed for IIR filters is poor since FILTERAUDIOFILE
%    must make use of the FILTER function instead of the FFT.
%  - For the FIR method, the function removes the delay introduced by 
%    the filter. That is not currently possible with the IIR method, since
%    I haven't figured out a way yet to use the FILTFILT function on 
%    consecutive windows while ensuring continuity between them.
%  - The IIR method is not recommended for time filtering of audio data
%    in chunks, with either overlap-save or overlap-add. Set ¦a¦ = 1 for
%    FIR filtering and calculate ¦b¦ using a function for FIR-filter
%    design (see FIR1 and FIRLS).
%
%  FUNCTION CALLS
%  1. xf = filterAudioFile(filePath,channel,b,a)
%     - For 'WAV' format. Variable input arguments FS, NCH, NBIT
%       and EN can also be included in the call but they'll be ignored
%       if FILEPATH has .wav extension.
%
%  2. xf = filterAudioFile(filePath,channel,b,a,nch,fs,nbit,en)
%     - Only for 'RAW' format
%
%  See also READAUDIOFILE, READWAVHEADER

%  VERSION 1.0
%  Guillermo Jimenez Arranz
%  email: gjarranz@gmail.com
%  11 Apr 2018

function xf = filterAudioFile(filePath,channel,samples,b,a,varargin)

[~,~,fext] = fileparts(filePath);

switch fext
    case {'.3gp','.aa','.aac','.aax','.act','.aiff','.amr','.ape',...
            '.au','.awb','.dct','.dss','.dvf','.flac','.gsm','.iklax',...
            '.ivs','.m4a','.m4b','.m4p','.mmf','.mp3','.mpc','.msv',...
            '.ogg','.oga','.mogg','.opus','.ra','.rm','.sln','.tta',...
            '.vox','.wma','.wv','.webm','.8svx'}
        error('Audio format not supported')

    case '.wav' % Import WAV file parameters

        % Error control: number of input arguments
        if nargin < 5, error('Not enough input arguments'); end
        if nargin > 5, error('Too many input arguments'); end

        % Read WAV header
        wavHeader = readwavHeader(filePath);

        % Size of audio data [bytes]
        fileStats = dir(filePath);
        if strcmp(wavHeader.subchunk2ID,'data')
            dataSize = wavHeader.subchunk2size;
            % subchunk2size = 0 when WAV file does not close properly (corrupt)
            if ~wavHeader.subchunk2size 
                dataSize = floor((fileStats.bytes-44)/wavHeader.blockAlign) ...
                    *wavHeader.blockAlign; % size of audio data [bytes]
                warning(['Size of audio data disagrees with header '...
                    'information: file is potentially corrupted. The '...
                    'data will be imported; note that this might result '...
                    'in a few initial samples being incorrectly read']);
            end
        elseif strcmp(wavHeader.subchunk2ID,'junk')
            dataSize = floor((fileStats.bytes-44)/wavHeader.blockAlign)...
                *wavHeader.blockAlign; % size of audio data [bytes]
        else
            error('Not recognised string for SUBCHUNK2ID')
        end

        % Audio parameters
        fs = wavHeader.sampleRate;
        nch = wavHeader.numChannels;
        nbit = wavHeader.bitsPerSample;
        en = wavHeader.byteOrder;

    otherwise % Import RAW file parameters
        % Error control: number of input arguments
        if nargin < 9, error('Not enough input arguments'); end
        if nargin > 9, error('Too many input arguments'); end

        % Size of audio data [bytes]
        filieStats = dir(filePath);
        dataSize = filieStats.bytes;

        % Audio parameters
        fs = varargin{1};
        nch = varargin{2};
        nbit = varargin{3};
        en = varargin{4};
end

% Determine the Section of Audio to be Read
    allSamples = dataSize*8/(nbit*nch);
    if ~isempty(samples)
        firstSample = min(max(samples(1),1),allSamples); % first audio sample to read (1 to ALLSAMPLES)
        nSamples = min(samples(2),allSamples - firstSample + 1); % number of audio samples to read (1 to ALLSAMPLES)
           
        % Error Control: Start Sample and Number of Samples
        if firstSample ~= samples(1)
            warning(['The starting sample SAMPLE(1) is out of bounds. '...
                'SAMPLE(1) = %d will be used'],firstSample)
        end        
        if nSamples ~= samples(2)
            warning(['The number of samples SAMPLE(2) is out of bounds. '...
                'SAMPLE(2) = %d will be used'],nSamples)
        end
    else
        firstSample = 1;
        nSamples = allSamples;
    end

% Zero-Pad Shorter Vector of Filter Coefficients
Na = length(a); % number of 'a' coefficients
Nb = length(b); % number of 'b' coefficients

% Assign Filter Type option
if Na > 1
    filtType = 'IIR'; 
else
    filtType = 'FIR'; 
end 
 
% Zero-Pad Shortest Vector of Coefficients (a,b)
N = max(Na,Nb);
a = [a(:) ; zeros(abs(N-Na),1)];
b = [b(:) ; zeros(abs(N-Nb),1)];

% Calculate Optimal Length of Overlapping Window (L)
nBitsSingle = 32; % number of bits for single-precision data (audio)
windowBytes = 52428800; % 50 MB window
nOverlap = N-1; % window overlap
Lt = nSamples + nOverlap;
L = 2^(floor(log2(windowBytes*8/nBitsSingle + nOverlap))); % window length (power of 2 for FFT efficiency)
if L > Lt, L = Lt; end % limit the maximum length of window to Lt
M = L - nOverlap; % length of audio chunk   

% Make M >= nOverlap to ensure only the first window is zero-padded
if M < nOverlap 
    error(['The filter order exceeds the maximum accepted for '...
        'efficient processing']);
end

% Filtering and Overlap-Save Parameters
nWindows = ceil(nSamples/M); % number of windows
zeropad = (nWindows)*M - nSamples; % length of zero-padding

% Initialise Variables
x_overlap = zeros(nOverlap,1); % overlap section
zi = zeros(1,nOverlap); % filter state       
s1 = 1;  % starting sample for the filtered vector xf

% FILTER SELECTED AUDIO
switch filtType
    % Overlap-Save Filtering (FILTER method)
    case 'IIR' 
        % Filter Audio Data
        hfig = waitbar(0,'','Name','Filtering'); % initialise waitbar
        for m = 1:nWindows-1
            % Read Audio Chunk
            iChunkSample1 = (m-1)*M + firstSample; % index of first sample to be read
            iChunkSample2 = m*M + firstSample - 1; % index of last sample to be read
            nChunkSamples = iChunkSample2 - iChunkSample1 + 1; % number of samples to read
            x_temp = readAudioFile(filePath,channel,...
                [iChunkSample1 nChunkSamples],'float',fs,nch,nbit,en); % audio chunk (without overlap)
            
            % Filter
            x = [x_overlap; x_temp]; % audio chunk (with overlap)
            x_overlap = x(M+1:L); % overlap-save for next audio chunk
            [xf_temp,zi] = filter(b,a,x,zi);
            clear x_temp x
            
            % Remove Initial Overlap (zero-pad) from First Window
            if m == 1
                xf_temp(1:round(nOverlap/2)) = [];
            end
            
            % Store Filtered Signal
            xf_temp = xf_temp(1:end - nOverlap);
            s2 = s1 + length(xf_temp) - 1;
            xf(s1:s2) = xf_temp;
            clear xf_temp
            
            % Parameters for Next Iteration
            s1 = s2 + 1;
            
            % Display processing status
            waitbarString = sprintf('Window %d/%d',m,nWindows);
            waitbar(m/nWindows,hfig,waitbarString); % show processing progress
        end

        % Read Audio Chunk
        iChunkSample1 = (nWindows-1)*M + firstSample;
        iChunkSample2 = firstSample + nSamples - 1;
        nChunkSamples = iChunkSample2 - iChunkSample1 + 1;
        x_temp = readAudioFile(filePath,channel,...
            [iChunkSample1 nChunkSamples],'float',fs,nch,nbit,en);
        
        % Filter
        x = [x_overlap ; x_temp ; zeros(zeropad,1)];
        xf_temp = filter(b,a,x,zi); % filter section padded with the last value
        clear x_temp x
        
        % Store Filtered Signal
        xf_temp = xf_temp(1:end - nOverlap - zeropad); % trim filtered chunk
        s2 = s1 + length(xf_temp) - 1;
        xf(s1:s2) = xf_temp;
        clear xf_temp
        
        % Display processing status
        waitbarString = sprintf('Window %d/%d',nWindows,nWindows);
        waitbar(1,hfig,waitbarString); % show processing progress
        delete(hfig) % delete waitbar

    % Overlap-Save Filtering (Fast-Fourier Transform Method)
    case 'FIR'    
        % Calculate the Fourier Transform of the Filter
        H = conj(fft(b,L)); % conjugate FFT of FIR filter b padded with M-1 zeros

        % Filter Audio Data
        hfig = waitbar(0,'','Name','Filtering'); % initialise waitbar
        for m = 1:nWindows-1
            % Read Audio Chunk
            iChunkSample1 = (m-1)*M + firstSample; % index of first sample to be read
            iChunkSample2 = m*M + firstSample - 1; % index of last sample to be read
            nChunkSamples = iChunkSample2 - iChunkSample1 + 1; % number of samples to read
            x_temp = readAudioFile(filePath,channel,...
                [iChunkSample1 nChunkSamples],'float',fs,nch,nbit,en);
            
            % Filter
            x = [x_overlap; x_temp];
            x_overlap = x(M+1:L);
            X = fft(x);
            xf_temp = ifft(X.*H);
            clear x_temp x X
            
            % Remove Initial Overlap (zero-pad) from First Window
            if m == 1
                xf_temp(1:round(nOverlap/2)) = [];
            end
            
            % Store Filtered Signal
            xf_temp = xf_temp(1:end - nOverlap); % trim filtered chunk
            s2 = s1 + length(xf_temp) - 1;
            xf(s1:s2) = xf_temp;
            clear xf_temp
            
            % Parameters for Next Iteration
            s1 = s2 + 1;
            
            % Display processing status
            waitbarString = sprintf('Window %d/%d',m,nWindows);
            waitbar(m/nWindows,hfig,waitbarString); % show processing progress
        end

        % Read Audio Chunk
        iChunkSample1 = (nWindows-1)*M + firstSample;
        iChunkSample2 = firstSample + nSamples - 1;
        nChunkSamples = iChunkSample2 - iChunkSample1 + 1;
        x_temp = readAudioFile(filePath,channel,...
            [iChunkSample1 nChunkSamples],'float',fs,nch,nbit,en);
        
        % Filter
        x = [x_overlap; x_temp; zeros(zeropad,1)];
        X = fft(x); % filter section padded with the last value
        xf_temp = ifft(X.*H);
        clear x_temp x X
        
        % Store Filtered Signal
        xf_temp = xf_temp(1:end - nOverlap - zeropad); % trim filtered chunk
        s2 = s1 + length(xf_temp) - 1;
        xf(s1:s2) = xf_temp;
        clear xf_temp
        
        % Display processing status
        waitbarString = sprintf('Window %d/%d',nWindows,nWindows);
        waitbar(1,hfig,waitbarString); % show processing progress
        delete(hfig) % delete waitbar
end
