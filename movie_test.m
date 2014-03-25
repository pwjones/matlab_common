
mr = VideoReader('/Users/pwjones/Movies/mouse_training/20121209/9085_2012-12-09-171006-0000.avi');
nFrames = 200;
vidHeight = mr.Height;
vidWidth = mr.Width;

% Preallocate movie structure.
mov(1:nFrames) = ...
    struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'),...
           'colormap', []);

% Read one frame at a time.

for k = 1 : nFrames
    mov(k).cdata = read(mr, k);
end

% Preallocate movie structure.
mov2(1:nFrames) = ...
    struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'),...
           'colormap', []);
       
for k = nFrame : (2*nFrames)
    mov2(k).cdata = read(mr, k);
end

% Size a figure based on the video's width and height.
hf = figure;
set(hf, 'position', [150 150 vidWidth vidHeight])

% Play back the movie once at the video's frame rate.
movie(hf, mov, 1, mr.FrameRate);