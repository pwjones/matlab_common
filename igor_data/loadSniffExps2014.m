% Load all of the processed data from the Spring-Fall 2014 sniffing dataset
save_flag = 1;
clear perMouseData;
base_folder = VIDEO_ROOT;
mouse_names = {'19439', '21413', '971', '1080', '1527'};
%mouse_names = {'19439'};
folders = {'140401', '140404', '140409', '141016', '140423', '140424', ...
           '141017', '141020', '141021', '141022', '141024', '141027', '141028', '141029'};
nMice = length(mouse_names);
following_thresh = 20; %mm
mov_thresh = 40; %mm/sec

allSniffData = [];
%for ii = 1:nMice
for ii = 3:3
    mouseData = struct('vids', [], 'resp', [], 'camTrig', []);
    for jj = 1:length(folders)
        saved_file = [base_folder filesep folders{jj} filesep folders{jj} '_' mouse_names{ii} '.mat'];
        if exist(saved_file, 'file')
            disp(['Loading: ' saved_file]);
            load(saved_file);
            
            mouseData.vids = cat(2, mouseData.vids, exp.vids);
            mouseData.resp = cat(2, mouseData.resp, exp.resp);
            mouseData.camTrig = cat(2, mouseData.camTrig, exp.camTrig);
        end
    end
    allSniffData = [allSniffData; mouseData];
    
    if save_flag
        save(mouse_names{ii}, 'mouseData', '-v7.3');
    end
end

if save_flag
    save('allSniffData', 'allSniffData', '-v7.3');
end
    
