% close all
% clear all

% select file
filedir=uigetdir();
cd(filedir);
filelist=dir('*.tiff');

%steping for temporal domain averaging
step=16;
start_frame=2400;
end_frame=4200;

tempimage1=zeros(480,640);

% determining k-space mask
% mask=zeros(640,640);
x=-320:319;
y=-320:319;
Y=y'*ones(1,length(x));
X=ones(length(y),1)*x;
kcenter=[0 98];
kradius=sqrt(kcenter(1).^2+kcenter(2).^2);

% reshape PSF

Ei=exp(1i.*2*pi.*(kcenter(1).*X/640+kcenter(2).*Y/640));
PSFre=zeros(640,640);
PSFre(321-191:321+219,321-162:321+118)=PSF;
shiftPSF=Ei.*PSFre;
shiftPSFfft=fftshift(fft2(shiftPSF,640,640));


rr=sqrt((X).^2+(Y).^2);

mask1=exp(-(rr-kradius).^2./10^2);
mask2=zeros(640,640);

%PSF masking

npixel=20;
for i=-319:320
    if i<kradius-npixel
        mask2(i+320,:)=1;
    else 
        mask2(i+320,:)=exp(-(i-npixel)^2/10^2);
    end
end

mask=mask1.*mask2;
maskPSF=mask.*shiftPSFfft;


%differential image

for j=0:step-1
        tempind=start_frame+j;
        temp_image=double(imread(filelist(tempind).name));
        tempimage1=tempimage1+temp_image;
end
tempimage1=tempimage1/step;

indlist=[];
intensity=[];
intensityper=[];
Eii=Ei(80:559,:);
frame_num=1;
new_im=zeros(480,640);


% iPM reconstruction

%calibration coefficient
coee=1e6;

for i=start_frame:step*1:end_frame
    
    tempimage2=zeros(480,640);
    
    for j=step:2*step-1
        tempind=i+j;
        temp_image=double(imread(filelist(tempind).name));
        tempimage2=tempimage2+temp_image;
    end
    tempimage2=tempimage2/step;
    
    diff_image=tempimage2-tempimage1;

    shift_diff_image=Eii.*diff_image;
    diff_fft=fftshift(fft2(shift_diff_image,640,640));

    mask_diff_fft=diff_fft.*mask;
    figure(3);imagesc(log10((abs(mask_diff_fft))));
    
    deconfft=mask_diff_fft./(maskPSF+1e-10);
    deconfft=deconfft.*mask;
    ss=find(abs(deconfft)>2);
    for kk=1:length(ss)
        ind=ss(kk);
        [rown,coln]=ind2sub([640 640],ind);
        deconfft(rown,coln)=(min(min(deconfft(max([1 rown-1]):min([640 rown+1]),max([1 coln-1]):min([640 coln+1]))))-deconfft(rown,coln))/((min([640 rown+1])-max([1 rown-1])+1)*(min([640 coln+1])-max([1 coln-1])+1)-1);
        if abs(deconfft(rown,coln))>1
            deconfft(rown,coln)=min(min(deconfft(rown-2:rown+2,coln-2:coln+2)));
        end

    end
    
    rec_diff=abs(fftshift(ifft2(deconfft)));
    rec_diff=rec_diff(1:480,:)*coee;

    figure(2);imagesc(rec_diff);colormap(hot);

    pause(0.05)
    frame_num=frame_num+1;
    
end
