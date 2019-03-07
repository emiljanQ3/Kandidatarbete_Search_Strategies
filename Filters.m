%% Fisheye correction for video
clc,clear

filepath='C:\Users\thoma\Desktop\Filmer\MTest';
listing=dir([filepath '*.mp4']);
for i = 1:2
FileName=listing(l).name;

V = VideoReader('011.mp4');
xsize = V.Height;
ysize = V.Width;
Frames = V.Duration*V.FrameRate;
mappingCoefficients = [880 -3e-4 0 0];       % mapping polynomial coefficients
imageSize = [xsize ysize];             % in [mrows ncols]
distortionCenter = [ysize./2 xsize./2];  % in pixels
intrinsics = fisheyeIntrinsics(mappingCoefficients,imageSize,distortionCenter);

video = VideoWriter('*m.avi');
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
%% Fisheye correction for picture
clc,clear

Barrel = imread('Barrel.jpg');
xsize = 1080;
ysize = 1920;
mappingCoefficients = [870 -3e-4 0 0];       % mapping polynomial coefficients
imageSize = [xsize ysize];             % in [mrows ncols]
distortionCenter = [ysize./2 xsize./2];  % in pixels
intrinsics = fisheyeIntrinsics(mappingCoefficients,imageSize,distortionCenter);
J = undistortFisheyeImage(Barrel,intrinsics);
%J = imcrop(J,[580,220,740,640]);
imshow(J)