clear all; close all;
% load the image file
[fname,mypath]=uigetfile({'*.tif;*.mat;*.fig'},'Select the image file (tif) ...');
stack=tiffread2([mypath fname]);
filename = find(mypath ==filesep);

filename = mypath((filename(end-1))+1:(filename(end))-1);
 


%extrack images from data structure
Total_Data = struct2cell(stack);
Images = Total_Data(7,1,:);
Images = cell2mat(Images);
L=size(Images);

% show projection of the stack and ROI
Image1(:,:) = max(Images,[],3);
f = figure;
colormap('default')
imagesc(Image1);
title('Select the cell body for analysis...double click');
h = imrect;
Neuron_position = wait(h);
title('Select section of background...double click');
h = imrect;
BG_position = wait(h);

close(f)
Sub_Im_Width = round(Neuron_position(3)/2);
Sub_Im_Length = round(Neuron_position(4)/2);
Sub_BG_Width = round(BG_position(3)/2);
Sub_BG_Length = round(BG_position(4)/2);

R_Neuron_X = uint16(Neuron_position(1) + Sub_Im_Width);
R_Neuron_Y = uint16(Neuron_position(2) + Sub_Im_Length);
R_Neuron_BG_X = uint16(BG_position(1) + Sub_BG_Width);
R_Neuron_BG_Y = uint16(BG_position(2) +Sub_BG_Length);

Neuron_Pixel_Size=100;
R_Neuron_Im = Images(R_Neuron_Y-Sub_Im_Length:R_Neuron_Y+Sub_Im_Length, R_Neuron_X-Sub_Im_Width:R_Neuron_X+Sub_Im_Width,:);
R_Neuron_Im = reshape(R_Neuron_Im,1,[],L(3));
R_Neuron_Im = squeeze(R_Neuron_Im);
R_Neuron_Im = R_Neuron_Im';
R_Neuron_BG_Im = Images(R_Neuron_BG_Y-Sub_BG_Length:R_Neuron_BG_Y+Sub_BG_Length, R_Neuron_BG_X-Sub_BG_Width:R_Neuron_BG_X+Sub_BG_Width,:);
R_Neuron_BG_Im = reshape(R_Neuron_BG_Im,1,[],L(3));
R_Neuron_BG_Im = squeeze(R_Neuron_BG_Im);
R_Neuron_BG = mean(R_Neuron_BG_Im);
[R_Neuron, IX]= sort(R_Neuron_Im, 2,'descend');
R_Neuron= R_Neuron(:,1:Neuron_Pixel_Size);
L=size(R_Neuron_Im);
R_Min_Col=R_Neuron(:,end);
R_Neuron_loc=zeros(L(1),L(2));

Right_Neuron_BG=zeros(L(1),L(2));

for i = 1:L(1)
   R_Neuron_loc(i,:)= R_Neuron_Im(i,:) >= R_Min_Col(i);
   Right_Neuron_BG(i,:)= R_Neuron_loc(i,:)*R_Neuron_BG(i);
end;

R_Neuron_Pixels=sum(R_Neuron_loc,2)';
Right_Neuron_BG=uint16(Right_Neuron_BG);
R_Neuron_loc=uint16(R_Neuron_loc);

R_Neuron_Im = R_Neuron_Im.*R_Neuron_loc;
CFP = R_Neuron_Im - Right_Neuron_BG;
CFP = sum(CFP,2)';
Norm_CFP=CFP./max(CFP);
Good_IND= find(Norm_CFP>=0.1);

CFP=CFP(Good_IND);

Right_Neuron_BG = sum(Right_Neuron_BG,2)'./R_Neuron_Pixels;
Right_Neuron_BG = Right_Neuron_BG(Good_IND);
 

Lab_View_File = load(strcat(mypath,fname(1:end-3),'txt'));
Temperat = Lab_View_File(:,2);
Temperature = interp1(1:length(Temperat),Temperat,Good_IND);
Image_Time = Good_IND;

baseline=mean(CFP(1:10));
delta_F=((CFP-baseline)/baseline)*100;

delta_F_mean = delta_F;
resampletemp = Temperature;
pathname = mypath;

Delta_R_Figure = figure;
axes_Delta_R = axes('Parent',Delta_R_Figure);
[AX,H1,H2]=plotyy(Image_Time,delta_F, Image_Time,Temperature,'plot');
hold on
%[AX,H1,H2]=plotyy(Image_Time,FRET1, Image_Time,Temperature,'plot');
set(H1,'LineStyle','-','LineWidth',2,'color', 'b');
set(H2,'LineStyle','--','LineWidth',2,'color', 'g');
%ylim(axes_Delta_R,[-50 150])
ylabel('\deltaF/F (%)');
xlabel('Time (sec) ');
ylabel(AX(2),'Temperature \circC')

saveas(Delta_R_Figure,[mypath filesep 'Delta_F_Figure_' filename '.fig'],'fig');
saveas(Delta_R_Figure,[mypath filesep 'Delta_F_Figure_' filename '.jpg'],'jpg');
filename = ['GCAMP_Analysis_' filename '.mat'];


        clear stack Total_Data Images Lab_View_File BG Image1 L_Neuron
                clear R_Neuron L_Neuron_BG_Im R_Neuron_BG_Im L_Neuron_Im
                clear L_Neuron_loc R_Neuron_Im R_Neuron_loc IX L_IX L_Neuron_BG
                clear L_Neuron_Pixels R_Neuron_BG R_Neuron_Pixels T_file M
                clear AX Delta_R_Figure H1 H2 L Neuron_Pixel_Size axes_Delta_R
                clear f

save ([mypath filename])

