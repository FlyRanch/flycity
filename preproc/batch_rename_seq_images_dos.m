% rename Photron image files (*.tif) if seqs mismatch
clear
clc
%warning off

dirs = dir;
for d=3:length(dirs)
    if dirs(d).isdir==1
        cd(dirs(d).name)
        subdirs=dir
        for s=3:length(subdirs)
            if subdirs(s).isdir==1
                cd(subdirs(s).name)




                    loc=cd;
%                     files = dir('*.bmp');
                    files = dir('*.tif');
                    infofile = dir('*.cih');
                    
                    if isempty(files)==0 && isempty(infofile)==0 

                        old_info = infofile(1).name;
                        new_info = [loc(end-12:end),old_info(end-3:end)];
                        if isequal(old_info,new_info)==0
%                             unix(['mv "' old_info '" "' new_info '"']);
                        %         movefile(old_info,new_info);
%                             java.io.File(old_info).renameTo(java.io.File(new_info));
                            dos(['rename "' old_info '" "' new_info '"']);
                        end

                        if isempty(files)==0
                            for i=1:length(files)

                                old_name = files(i).name;
                                new_name = [loc(end-12:end),old_name(end-9:end)];

                                if isequal(old_name,new_name)==0
%                                     unix(['mv "' old_name '" "' new_name '"']);
                            %         copyfile(old_name,new_name);
                        %             movefile(old_name,new_name);
                            %         java.io.File([old_name]).renameTo(java.io.File([new_name]));
                                    dos(['rename "' old_name '" "' new_name '"']);
                                end

                                counter=length(files)-i
                            end
                        end
                    end
                    
                    cd ..
                    
            end
        end
        cd ..
    end
end

    