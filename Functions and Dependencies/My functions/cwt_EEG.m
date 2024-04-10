function [EEGTensor, sigma] = cwt_EEG(data, fs)
%CWT_EEG creates an EEG tensor from the input data using the CWT.
%   This function applies the CWT using a mexican hat wavelet aplied to the
%   EEG data. The EEG data is first normalized before applying the CWT.

%
%INPUT
%   data : a (channels x time) matrix
%   fs   : the sampling frequency
%
%OUTPUT
%   EEGTensor : a (channels x time x scales) tensor
%   sigma     : the std deviation of the original data 

%This function is writen by the course staff of the course EE4750 and is only slightly altered
%by Mees van Dijk.

%% Amount of scales can be adjusted here, the default is 30 scales.
steps = 30;
scale_steps = linspace(1,30,steps);
%% normalize the data per channel
m = mean(data,2); 
data = data-m(:,ones(1,size(data,2)));

sigma = std(data,[],2);
data = data./sigma(:,ones(1,size(data,2)));

%% Apply wavelet transform on the data per channel
wl_name = 'mexh';                           % The type of wavelet we use - here Mexican hat
wl_scales =  make_scales(scale_steps,wl_name,fs);  % The scales at which the cwt is calculated (here: compute scales for frequencies between 1 and 30 Hz)

EEGTensor = zeros(size(data,1),steps,size(data,2));
for channel = 1:size(data,1)
    cwtdata = cwtft(data(channel,:),'scales',wl_scales,'wavelet',wl_name);
    EEGTensor(channel,:,:) = cwtdata.cfs;          % channel-time-frequency data
end

end

function scales = make_scales(freq_vector,wl,fs) 
% scales = make_scales(freq_vector,wl,fs) 
% To create the wavelet scales corresponding to the frequencies determined
% by freq_vec
%% INPUTS:
% freq_vector: the vector of frequencies
% wl: the name of the wavelet function
% fs: sampling frequency
%% OUTPUTS:
% scales: wavelet scales

scales=zeros(1,numel(freq_vector));
centerFreq = centfrq(wl);
for i=1:numel(freq_vector)
    scales(i) = centerFreq*fs/freq_vector(i);
end
end