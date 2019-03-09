%% Filter for videofile
clc,clear
tic

filepath = 'C:\Users\thoma\Desktop\Filmer\KKV\';
listing = dir([filepath '*.mp4']);
for i = 1:1:length(listing)
VideoName = listing(i).name;

V = VideoReader(VideoName);
xsize = V.Height;
ysize = V.Width;
Frames = V.Duration*V.FrameRate;
mappingCoefficients = [880 -3e-4 0 0];       % mapping polynomial coefficients
imageSize = [xsize ysize];             % in [mrows ncols]
distortionCenter = [ysize./2 xsize./2];  % in pixels
intrinsics = fisheyeIntrinsics(mappingCoefficients,imageSize,distortionCenter);

video = VideoWriter(VideoName);
video.FrameRate = V.Framerate;
open(video);

while hasFrame(V)
    FrameN = readFrame(V);
    CorrectedFrame = undistortFisheyeImage(FrameN,intrinsics);  % Remove fisheye
    CorrectedFrame = imcomplement(CorrectedFrame);  % Invert
    CorrectedFrame = imcrop(CorrectedFrame,[580,220,740,640]);  % Cropping
%    CorrectedFrame = rgb2gray(CorrectedFrame);  % Grayscale
    writeVideo(video,CorrectedFrame);
end

close(video)
end
disp('Dun')
toc
%% Filter for single video
clc,clear
tic

V = VideoReader('011.mp4');
xsize = V.Height;
ysize = V.Width;
Frames = V.Duration*V.FrameRate;
mappingCoefficients = [880 -3e-4 0 0];       % mapping polynomial coefficients
imageSize = [xsize ysize];             % in [mrows ncols]
distortionCenter = [ysize./2 xsize./2];  % in pixels
intrinsics = fisheyeIntrinsics(mappingCoefficients,imageSize,distortionCenter);
video = VideoWriter('011m.avi');
video.FrameRate = V.Framerate;
open(video);

while hasFrame(V)
    FrameN = readFrame(V);
    CorrectedFrame = undistortFisheyeImage(FrameN,intrinsics);  % Remove fisheye
    CorrectedFrame = imcomplement(CorrectedFrame);  % Invert
    CorrectedFrame = imcrop(CorrectedFrame,[580,220,740,640]);  % Cropping
%    CorrectedFrame = rgb2gray(CorrectedFrame);  % Grayscale
    writeVideo(video,CorrectedFrame);
end

close(video)

disp('Dun')
toc
%% Fisheye correction for picture
clc,clear

Barrel = imread('Barrel.jpg');
xsize = 1080;
ysize = 1920;
mappingCoefficients = [880 -3e-4 0 0];       % mapping polynomial coefficients
imageSize = [xsize ysize];             % in pixels
distortionCenter = [ysize./2 xsize./2];  % in pixels
intrinsics = fisheyeIntrinsics(mappingCoefficients,imageSize,distortionCenter);
J = undistortFisheyeImage(Barrel,intrinsics);
%J = imcrop(J,[580,220,740,640]);
imshow(J);