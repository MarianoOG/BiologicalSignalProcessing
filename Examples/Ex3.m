% Activity 3: design of a digital filter from an analog filter
% Design a band-pass filter from 20 Hz to 500 Hz using a low pass and a 
% high pass in cascade
% 1) Obtain H(s)
% 2) Obtain H(z) using bilinear transormation
% 3) Sketch H(z) for each filter use freqz command.
figure
freqz([0.6498,0.6498],[2.6498,-1.3502],[])
figure
freqz([2,-2],[2.0251,-1.9749],[])
figure
freqz([1.2996,0,-1.2996],[5.3661,-7.9674,2.6665],[])