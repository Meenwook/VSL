%% VLC Receiver Code, Author: Meenwook Ha

clc;clear;

imgTime=0.3; %Input the image-time for transmitter (The time an image is lasted).

endpoint=100; %Set the number of frames I will take  

c1=zeros(7,endpoint); c2=zeros(7,endpoint); %Initialize 2 Hamming matrices

P=[1 1 1 0;0 1 1 1;1 1 0 1]; 
H=horzcat(P,eye(3)); %Parity Check Matrix

%%1. Basic Setting for webcam.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vid=imaq.VideoDevice('winvideo',1);
set(vid.DeviceProperties,'HorizontalFlip','on');
set(vid.DeviceProperties,'BacklightCompensation','off');
set(vid.DeviceProperties,'ExposureMode','manual','FocusMode','manual');
set(vid.DeviceProperties,'Brightness',100,'Gain',60,'WhiteBalance',2000);
set(vid.DeviceProperties,'Focus',0,'FrameRate','10.0000');
frameRate=vid.DeviceProperties.FrameRate; %10 FPS is the max FPS
frameRate=str2double(frameRate);
samplingR=frameRate*imgTime; %Calculate sampling rate for analysis

%Set the ranges of red and blue blob sizes
blobB = vision.BlobAnalysis('AreaOutputPort', false,...
'MinimumBlobArea', 300,'MaximumBlobArea', 100000,'MaximumCount',6);
blobR = vision.BlobAnalysis('AreaOutputPort', false,...
'MinimumBlobArea', 150,'MaximumBlobArea', 100000,'MaximumCount',14);

%Set box shaper so that it is possible to know if the frames are well
%recognized while communication.
boxShpB = vision.ShapeInserter('Fill',true,'FillColorSource','Input port','Opacity',0.8);
boxShpR = vision.ShapeInserter('Fill',true,'FillColorSource','Input port','Opacity',1);
vidP = vision.VideoPlayer('Name', 'Detector','Position', [100 100 640 480]);



%%2. Code for Communication%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic; %Count Communication Time
for frameNum=1:1:endpoint

rbF = uint8(255*step(vid)); %Take a snapshot 10FPS

[fbinR,fbinB]=maskRB(rbF); %Detect Colors using color threshold function
[lcB,xB]=step(blobB,fbinB);lcB=uint16(lcB); %Insert Blue first

detF=rbF;

%Process below detection process below iff all boundaries (6) are on.
if length(lcB(:,1))==6 && length(lcB(:,1))>=1
avgBH=mean(lcB(:,2));

[lcR,xR]=step(blobR,fbinR);lcR=uint16(lcR);

%Detect Red Digits and their coordinates.
for k=1:length(lcR(:,1))
    if lcR(k,1)<lcB(1,1)
        if lcR(k,2)>avgBH;c1(7,frameNum)=1;
        else;c2(7,frameNum)=1;
        end
    elseif lcR(k,1)>lcB(1,1) && lcR(k,1)<lcB(2,1)
        if lcR(k,2)>avgBH;c1(6,frameNum)=1;
        else;c2(6,frameNum)=1;
        end
    elseif lcR(k,1)>lcB(2,1) && lcR(k,1)<lcB(3,1)
        if lcR(k,2)>avgBH;c1(5,frameNum)=1;
        else;c2(5,frameNum)=1;
        end
    elseif lcR(k,1)>lcB(3,1) && lcR(k,1)<lcB(4,1)
        if lcR(k,2)>avgBH;c1(4,frameNum)=1;
        else;c2(4,frameNum)=1;
        end
    elseif lcR(k,1)>lcB(4,1) && lcR(k,1)<lcB(5,1)
        if lcR(k,2)>avgBH;c1(3,frameNum)=1;
        else;c2(3,frameNum)=1;
        end
    elseif lcR(k,1)>lcB(5,1) && lcR(k,1)<lcB(6,1)   
        if lcR(k,2)>avgBH;c1(2,frameNum)=1;
        else;c2(2,frameNum)=1;
        end
    else 
        if lcR(k,2)>avgBH;c1(1,frameNum)=1;
        else;c2(1,frameNum)=1;
        end
    end    
end

% Insert BoxShaper in the output screen
detF=step(boxShpR,detF,xR,uint8([255 0 0]));
detF=step(boxShpB,detF,xB,uint8([0 0 255]));
end
    detF = im2single(detF);
    step(vidP, detF)
end

timeElapsed = toc;

%%3. Analyzing Part%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Error Detection using Syndrome vectors
s1=mod(H*c1,2);s2=mod(H*c2,2);
c1Synd=c1; c2Synd=c2;

for k=1:length(s1(1,:))
    for kk=1:7
    if isequal(s1(:,k),H(:,kk))==1
       if c1Synd(kk,k)==1; c1Synd(kk,k)=0;
       else; c1Synd(kk,k)=1;
       end
    end
    end
end

for k=1:length(s2(1,:))
    for kk=1:7
    if isequal(s2(:,k),H(:,kk))==1
       if c2Synd(kk,k)==1;c2Synd(kk,k)=0;
       else;c2Synd(kk,k)=1;
       end
    end
    end
end

%Except parity bits and connect two codebooks
cSynd=vertcat(c1Synd(1:4,:),c2Synd(1:4,:));

%Find Pilot Vector and truncate based on pilot vectors
pV=[1;1;1;1;1;1;1;1];
pDV=zeros(1,length(cSynd));

for k=1:length(cSynd)
    if isequal(cSynd(:,k),pV)
    pDV(1,k)=1;
    end
    
end

pL=find(pDV==1);
pDif=diff(pL);
[pDifV]=max(pDif);
firstPEnd=find(pDif==pDifV);
startP=pL(firstPEnd)+1;
endP=pL(firstPEnd+1)-1;
cTrunc=cSynd(:,startP:1:endP);

cTrans=cTrunc';

%%Synchronization Error detection loop
cSync=cTrans;
[~,~,z]=unique(cSync,'row','legacy');
befoSync=find(diff(z));
aftSync=befoSync;
k=1;

while k<=length(aftSync)
    if k>length(aftSync);break;end
    if rem(aftSync(k),samplingR)==0         %If no error, go to next vector
       k=k+1;   
    elseif rem(aftSync(k),samplingR)==1     %If an error or oversampled
                                            %, remove it
            cSync(aftSync(k),:)=[];
    elseif rem(aftSync(k),samplingR)==2     %If undersampled, add one vector
            cSync=[cSync(1:aftSync(k),:)
                   cSync(aftSync(k),:)
                   cSync(aftSync(k)+1:end,:)];

        
    end
    [~,~,z]=unique(cSync,'row','legacy');
            aftSync=find(diff(z));   
end

%Detection for the last information vector
cEdited=cSync;
    if rem(length(cEdited),samplingR)==1
        cEdited(length(cEdited(:,1)),:)=[];
    elseif rem(length(cEdited),samplingR)==2
        cEdited=[cEdited(1:length(cEdited),:)
                   cEdited(length(cEdited),:)];
    end
        
%Downsample & Output
alpLength=length(cEdited)/samplingR;
cDown=zeros(alpLength,8);
for k=1:alpLength
cDown(k,:)=cEdited(k*samplingR,:);
end
cFlip=fliplr(cDown);
cBin=bi2de(cFlip);
cOut=zeros(1,alpLength);
for k=1:length(cBin)
    cOut(k)=char(cBin(k));
end
fprintf('\n The Received word is: %s\n',cOut);

