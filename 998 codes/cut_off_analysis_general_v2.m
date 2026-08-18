function [fourier_position,fourier_intensity,dark_position]=cut_off_analysis_general(main_folder,start_pos,step_pos,figShow)


intensity=[];

img_sv=[];

%main_folder='F:\2024\20240521_BMOL\Data\005_8340eV\002_cutoff_scan\';



x=[];

%start_pos=0.750267;
%step_pos=-0.001;
kk=1;

original_files=dir([main_folder '/*.tif']); 
dark_position=[main_folder '/' original_files(2).name];
dark=0.*imread(dark_position);

 for kk=2:length(original_files)

    img=imread([main_folder '/' original_files(kk).name]);
    tokens = regexp(original_files(kk).name, '[-+]?\d+\.?\d*', 'match');
    cutoff_pos = str2num(tokens{3});
    
    img=img-dark;
    %figure; imagesc(img)
    %img=circshift(img,[200 0]);
    %img=mcrop(img, [500,500]);
    img=mcrop(img, min(size(img),[500 500]));
    
    intensity(end+1)= mean(img(:));
    %size(img_sv)
    img_sv(:,:,end+1)=img;
    
    %size(cutoff_pos)
    x(end+1)=cutoff_pos;%start_pos+(kk-1)*step_pos;
    %kk=kk+1;
end
%figure; imagesc(img); axis image; colormap gray;
%figure; plot(x,intensity);
%figure; plot(x(1:end-1),-diff(intensity));


%x

[uniqueA unique_i j] = unique(x,'first');
%indexToDupes = find(not(ismember(1:numel(x),i)))
x=x(unique_i);
intensity=intensity(unique_i);

%intensity


xq=x(1):(x(end)-x(1))/(length(original_files)-1):x(end);
x(1)
x(end)
(x(end)-x(1))/(length(original_files)-1)

intensity2 = interp1(x(:),intensity(:),xq(:));
if figShow
    figure; plot(x,intensity); hold on;plot(xq,intensity2);% hold on; plot(intensity2);
end
%intensity2

%fourier_position=x(1:end-1);
%fourier_intensity=diff(intensity);
fourier_position=xq(1:end-1);
fourier_intensity=diff(intensity2);

%figure; plot(fourier_intensity);
%figure; plot(intensity);

%fourier_position=flip(fourier_position);
%fourier_intensity=flip(fourier_intensity);

%% gaussian fitting to extend the edges 
%{
add=8;
fourier_position_extended=fourier_position;
fourier_position_extended=cat(2,(-add:-1).*(fourier_position_extended(2)-fourier_position_extended(1))+fourier_position_extended(1),fourier_position_extended);
fourier_position_extended=cat(2,fourier_position_extended,(1:add).*(fourier_position_extended(2)-fourier_position_extended(1))+fourier_position_extended(end));
%data_extended=-exp(-29000*(fourier_position_extended-5.5775).^2+5.6);

% design matrix for least squares fit
xdata = fourier_position(:);
A = [xdata.^2,  xdata,  ones(size(xdata))]; 

% log of your data 
b = log(-fourier_intensity(:));                  

% least-squares solution for x
x = A\b;

data_extended=-exp(x(1)*(fourier_position_extended).^2+x(2)*(fourier_position_extended)+x(3));

figure; plot(fourier_position,(fourier_intensity)); hold on; plot(fourier_position_extended,data_extended);

fourier_position=fourier_position_extended;
fourier_intensity=data_extended;
%}
%% low pass
%{
%low_pass=1.8/5;
low_pass=5;
%low_pass=1/4;
fourier_intensity=cat(1,fourier_intensity(:),flip(+fourier_intensity(:),1));
%figure; imagesc(fourier_intensity)
fourier_intensity=s_fft(fourier_intensity,1);

kkx=(1:size(fourier_intensity,1))-floor(size(fourier_intensity,1)/2)-1;kkx=kkx./size(fourier_intensity,1);
kkx=reshape(kkx,[],1);
filter=exp(-(low_pass*kkx).^2);

fourier_intensity=fourier_intensity.*filter;
%figure; plot(abs(fourier_intensity));

%fourier_intensity(1:round(end*low_pass))=0;fourier_intensity(end-round(end*low_pass):end)=0;
fourier_intensity=real(s_ifft(fourier_intensity,1));
fourier_intensity=fourier_intensity(1:end/2);
%fourier_intensity=fourier_intensity-mean(fourier_intensity(end-5:end));
%fourier_intensity(fourier_intensity>0)=0;
%}
%% gaussian fitting

%fourier_position=fourier_position(20:60);
%fourier_intensity=fourier_intensity(20:60);

%fourier_intensity(1)=30.*fourier_intensity(1)

fourier_intensity=+fourier_intensity;

%figure; plot(fourier_position,fourier_intensity);
%figure; plot(x,intensity);



%save([main_folder 'result'],'fourier_position','fourier_intensity')

