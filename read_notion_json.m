function [data, t_full, deviceNames] = read_notion_json(expType,expDate)
exp_day_str = datestr(expDate,'yyyymmdd');
ext = '.json';
baseDir = 'C:\Users\BatLab\Documents\eeg_data\';
expDir = dir(fullfile(baseDir,['eeg_recording_' exp_day_str '*']));
deviceNames = {'776','953'};
nDevice = length(deviceNames);
n_samp_per_chunk = 25;
fs = 250;
nFile = length(dir(fullfile(expDir.folder,expDir.name,[expType '*' ext])))/length(deviceNames);
[data,t_full] = deal(cell(nFile,nDevice));

for d_k = 1:nDevice
    fNames = dir(fullfile(expDir.folder,expDir.name,[expType '*' deviceNames{d_k} '*' ext]));
    t = deal(cell(nFile,1));
    for f_k = 1:nFile
        str = fileread(fullfile(fNames(f_k).folder,fNames(f_k).name));
        idx = strfind(str,'data');
        idx = [idx length(str)+3];
        for k = 1:length(idx)-1
            val(k) = jsondecode(str(idx(k)-2:idx(k+1)-3));
        end
        info = [val.info];
        t{f_k} = [info.startTime];
        data{f_k,d_k} = cat(2,val.data)';
    end
    
    for f_k = 1:nFile
        t_interp = nan(length(t{f_k}),n_samp_per_chunk);
        for chunk_k = 1:length(t{f_k})
            t_interp(chunk_k,:) = t{f_k}(chunk_k) + linspace(0,1e3*(n_samp_per_chunk-1)/fs,n_samp_per_chunk);
        end
        t_full{f_k,d_k} = 1e-3*reshape(t_interp',1,[])';
    end
end
end