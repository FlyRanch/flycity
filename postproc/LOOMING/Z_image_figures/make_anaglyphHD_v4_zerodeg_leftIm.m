% make anaglyph from L&R
clear
clc

paralax = 1;

cutx = 800
cuty = 500
rez=300; %resolution (dpi) of final graphic

list_left = dir('3DrenderLeft*fig');
% list_right = dir('3DrenderRight*fig');

for i = 1:length(list_left)
    nr = list_left(i).name(13:18);
    
    open(list_left(i).name)
    f=gcf;
    
    view(0-paralax,0)
%     saveas(gca,['3DrenderLeft.tif'])
    print(f,'3DrenderLeft.tif','-dpng',['-r',num2str(rez)],'-opengl') %save file 
%     print(f,'3DrenderLeft.tif','-dtiffn',['-r',num2str(rez)],'-opengl') %save file 

    view(0+paralax,0)
%     saveas(gca,['3DrenderRight.tif'])
    print(f,'3DrenderRight.tif','-dpng',['-r',num2str(rez)],'-opengl') %save file 
%     print(f,'3DrenderRight.tif','-dtiffn',['-r',num2str(rez)],'-opengl') %save file 
    close(f)
    
%     open(list_right(i).name)
%     view(0+paralax,0)
%     saveas(gca,['3DrenderRight.tif'])
%     f=gcf;
% %     print(f,'3DrenderRight.tif','-dtiffn',['-r',num2str(rez)],'-opengl') %save file 
%     close(f)

 leftEyeImage  = imread(['3DrenderLeft.tif'])  ; % Load the left eye image.
 rightEyeImage = imread(['3DrenderRight.tif']) ; % Load the right eye image.
 
 leftEyeImageCrop = leftEyeImage(cuty:end-cuty,cutx:end-cutx,:);
 imwrite(leftEyeImageCrop,['3DrenderLeft_crop_0deg_fr' nr '.tif'])
 
 leftEyeImage(:,:,2:3) = 0 ;               % Removes green and blue from the left eye image.
 rightEyeImage(:,:,1)  = 0 ;               % Removes red from the right eye image.
 anaglyph = leftEyeImage + rightEyeImage ; % Combines the two to produce the finished anaglyph.
 
 anaglyph = anaglyph(cuty:end-cuty,cutx:end-cutx,:);
 
 imwrite(anaglyph,['3DrenderAnaglyph_crop_0deg_fr' nr '.tif'])

%  
%  imshow(anaglyph(cut+1:end-cut,2*cut+1:end-2*cut,:),'border','tight') ;       % Show the anaglyph image with no padding.
% %  print(gcf,'-dtiffn','-painters',['3DrenderAnaglyphHD' nr '.tif'])  % Save the anaglyph image.
%     f=gcf;
%     print(f,['3DrenderAnaglyph_fr' nr '.tif'],'-dtiffn','-painters',['-r',num2str(rez)],'-opengl') %save file 
%     close(f)
end

