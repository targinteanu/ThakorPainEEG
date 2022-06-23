%% MODULE 2 
%% SECTION 1: Bandpass filter the EEG data

% extract data from EEG structure
data = EEG.data;

% TODO: Fill in fs - the sampling frequency.
% HINT: EEG is the structure that stores all recorded EEG information.
fs = EEG.srate; % the signal sampling frequency (in Hz)

numChan = size(data,1);
fData = zeros(size(data));

% TODO: Fill in lowpass and highpass frequencies
fl = .1; % the lower bound cutoff frequency (in Hz)
fh = 170;% the higher bound cutoff frequency (in Hz)

% bandpass filter each channel of the recorded EEG data
for i = 1:numChan
    fData(i,:) = bpf(double(data(i,:)), fl, fh, fs);
end
EEG.data = single(fData);

%% SECTION 2: Plot the spectral content of first EEG channel before and after bandpass filtering
figure();
titles = { 'RAW EEG', sprintf( 'BPF EEG (%.2f - %.2f Hz)', fl, fh ) };
channels = [ data(1,:); EEG.data(1,:) ];
for i = 1:2
    signal = channels(i,:);
    signal_length = length( signal );           % calculate the signal length
    signal = detrend( signal );                 % remove the linear trend of the signal
    fft_length = 2 ^ nextpow2( signal_length ); % compute the next power of 2 from the signal length
    
    % TODO: Fill in the code to compute the FFT and frequency range
    % HINT: Use the in-built 'fft' function
    Y = fft(signal);% the FFT of the EEG channel
    f = linspace(0, fs/2, fft_length/2); % the frequency range of the single-sided FFT
    
    % plot the single-sided FFT
    subplot( 2, 1, i )
    plot( f, 2 * abs( Y(1:fft_length/2) ) );
    
    xlabel( 'Frequency (Hz)' );
    ylabel( '|Y(f)|' );
    title( titles{ i } );
end

%% SECTION 3: Save the filtered .set file
pop_saveset(EEG);

%% SECTION 4: Generate topological maps of EEG activity
EEG = pop_loadset( char( 'eeg_ica.set' ) );     % load the ICA dataset

fft_winlen = 0.25;      % window length for the FFT

% TODO: Determine the frequency band ranges (in Hz) for the Delta, Theta, Alpha, Beta, and Gamma bands
% NOTE: Use the following format: [ delta_low   delta_high;
%                                   theta_low   theta_high;
%                                   alpha_low   alpha_high;
%                                   beta_low    beta_high;
%                                   gamma_low   gamma_high ]
freq_band = [0.5, 4;
             4, 7;
             8, 12;
             12.5, 30;
             25, 140];

% extract a window of time 450 ms to 850 ms after the stimulus was presented
eeg_window = double( EEG.data(1:62, 725:925, :) );

for i = 1:size( eeg_window, 3 ) % for each trial (10 in total)
    data = eeg_window( :, :, i );
    [ spectra, F, T ] = getspectralpower( data, fft_winlen * fs, 0, fft_winlen * fs, fs );
    [ nc, nf, nt ] = size( spectra );
    
    bandnum = size( freq_band, 1 );	% number of frequency bands we are computing
    bandpower = nan( nc, bandnum ); % preallocate the signal power in each band for each channel
    
    for c = 1:nc            % for every channel
        for b = 1:bandnum   % for every frequency band
            bandpower( c, b ) = squeeze( mean( spectra( c, F >= freq_band( b, 1 ) & F <= freq_band( b, 2 ) ), 2 ) );
        end
    end
    
    % compute the relative power in each frequency band
    relative_power = zeros( nc, bandnum );  % preallocate the relative power for each frequency band
    bandpower_sum = sum( bandpower, 2 );    % compute the power sum for each frequency band
    for c = 1:nc            % for each channel
        for b = 1:bandnum   % for every frequency band
            relative_power( c, b ) = bandpower( c, b ) ./ bandpower_sum( c, 1 );
        end
    end
    bandpower_total( :, :, i ) = relative_power;
end

% plot the topographical map for each power band
titles = { sprintf( 'Delta (%.2f - %.2f Hz)', freq_band( 1, 1 ), freq_band( 1, 2 ) ), ...
           sprintf( 'Theta (%.2f - %.2f Hz)', freq_band( 2, 1 ), freq_band( 2, 2 ) ), ...
           sprintf( 'Alpha (%.2f - %.2f Hz)', freq_band( 3, 1 ), freq_band( 3, 2 ) ), ...
           sprintf( 'Beta (%.2f - %.2f Hz)',  freq_band( 4, 1 ), freq_band( 4, 2 ) ), ...
           sprintf( 'Gamma (%.2f - %.2f Hz)', freq_band( 5, 1 ), freq_band( 5, 2 ) ) };
bandpower_avg = mean( bandpower_total, 3 ); % take the average over all trials

figure();
for i = 1:bandnum
    subplot( 1, bandnum, i );
    topoplot( bandpower_avg( :, i ), EEG.chanlocs, 'maplimits', 'maxmin' );
    title( titles{ i } );
    colorbar;
end