function [fourier_position,fourier_intensity,dark_position]=cut_off_analysis_general(main_folder,start_pos,step_pos)


intensity=[];

img_sv=[];

%main_folder='F:\2024\20240521_BMOL\Data\005_8340eV\002_cutoff_scan\';



x=[];

%start_pos=0.750267;
%step_pos=-0.001;
kk=1;

original_files=dir([main_folder '/*.tif']); 
dark_position=[main_folder '/' original_files(1).name];
dark=imread(dark_position);

 for kk=1:length(original_files)

    img=imread([main_folder '/' original_files(kk).name]);
    tokens = regexp(original_files(kk).name, '(\d+)\.(\d+)', 'tokens', 'once');
    numbers = str2double(vertcat(tokens{:}));
    numbers
    img=img-dark;
    %figure; imagesc(img)
    %img=circshift(img,[200 0]);
    img=mcrop(img, [800,800]);
    intensity(kk)= mean(img(:));
    %size(img_sv)
    img_sv(:,:,kk)=img;
    x(end+1)=start_pos+(kk-1)*step_pos;
    kk=kk+1;
end
%figure; imagesc(img); axis image; colormap gray;
%figure; plot(x,intensity);
%figure; plot(x(1:end-1),-diff(intensity));

fourier_position=x(1:end-1);
fourier_intensity=diff(intensity);


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
low_pass=1.8/5;
%low_pass=1/4;
fourier_intensity=cat(1,fourier_intensity(:),flip(+fourier_intensity(:),1));
%figure; imagesc(fourier_intensity)
fourier_intensity=s_fft(fourier_intensity,1);

fourier_intensity(1:round(end*low_pass))=0;fourier_intensity(end-round(end*low_pass):end)=0;
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

