function Background = GetBackgroundv1(frame,PAR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[Program Name]
%  GetBackgroundv1.m
%[Description]
%  Generate background images
%[Usage]
%  GetBackground( exp002,frame)
%[Programmer]
%  Ebraheem Fontaine
%[Development Environment]
%  Matlab R2007a
%[Version]
%  ver.1.0  06/08/2007 Initial Version (exp002)
%  ver.1.1	06/12/2007 exp002, exp035, exp83, exp098, exp101
%  ver.1.0  01/14/2008 Allow user to select the region to smooth over and
%  make the background
%[Input]
%  dir_name     Input and output directory name
%  exp_name     Experiment name
%[Return]
%  Background   background image
%[Output]
%  Background images (jpeg files)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[variable Name]
%  cam                      Camera number
%  frame                    Frame number
%  p_frame                  Frame number for iteration
%  input_filename           Input file name
%  output_filename          Output file name
%  Input(i,j,cam)           Original (input) image
%  Image(i,j,cam)           Image after contrast enhancement
%  Background(i,j,cam)      Background image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for cam = 1:3

    % Load images for background image generation
    % Camera 1 (Frame 1)
    
%     %Number of digit places for the total number of frames
%     % photron file structure
%     digits = 6;
% %     digits = length(num2str(PAR.numframes));
%     
%     input_filename = sprintf(['%scam%03d/%s%0' num2str(digits) 'd%s'], ...
%         PAR.imagepath, cam, PAR.stub, frame,PAR.image_filter(2:end));
%     

    % photron file struc
    input_filename = ([PAR.imagepath, PAR.stub(1:3), num2str(cam), PAR.stub(5:end), '/', PAR.stub(1:3), num2str(cam), PAR.stub(5:end), sprintf(['%0',num2str(3),'d'], frame),PAR.image_filter(2:end)]);
    
    Input = imread(input_filename);

%% no BG contrast enhancement yet, but in filter file 'FlySegment.m' (FTMmod 20120605)
%     % Contrast enhancement
%     Image = imadjust(Input);
%     
%     imagesc(Image); colormap gray
    
    imagesc(Input); colormap gray
    
    % remove fly?
    rem_fly = input('filter fly from image? Y/N [N]: ', 's');
    if isempty(rem_fly) == 1
%         Background(:,:,cam) = Image;
        Background(:,:,cam) = Input;
    elseif rem_fly == 89 || rem_fly == 121
        Background(:,:,cam) = roifill;
    else
%         Background(:,:,cam) = Image;
        Background(:,:,cam) = Input;
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
