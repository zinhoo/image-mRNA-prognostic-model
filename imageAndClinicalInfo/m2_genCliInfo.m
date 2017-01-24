% Generate clinical info. Patients with survival time not available or
% less than 30 days, or without vital status, are removed. Patients whose
% grade and stage are not available are also removed

% cliInfo.pid:
% cliInfo.age;
% cliInfo.gender: 'MALE' or 'FEMALE'
% cliInfo.grade: 'G1', 'G2', 'G3', or 'G4'
% cliInfo.stage: 'State I', 'Stage II', 'Stage III', or 'Stage IV'
% cliInfo.time
% cliInfo.death

clear

fid = fopen('nationwidechildrens.org_clinical_patient_kirc.txt');


% Discard the first three rows.
for i = 1 : 3
    line = fgetl(fid);
end

n = 0;
line = fgetl(fid);
while ischar(line)
    if isempty(strtrim(line))
        break;
    end   
    
    c = regexp(line, '\t', 'split');
    c2 = regexp(c{2}, '-', 'split');
    
    if (strcmp(c{25}, 'Alive') || strcmp(c{25}, 'Dead'))
        n = n + 1;
        cliInfo.pid{n, 1} = [c{2}, '-01']; 
        cliInfo.age(n, 1) = str2double(c{45}); 
        cliInfo.gender{n, 1} = c{9};   
        cliInfo.grade{n, 1} = c{5}; 
        cliInfo.stage{n, 1} = c{24}; 
        if strcmp(c{25}, 'Alive')
            cliInfo.time(n, 1) = str2double(c{26});
            cliInfo.death(n, 1) = 0;
        else
            cliInfo.time(n, 1) = str2double(c{27});
            cliInfo.death(n, 1) = 1;
        end            
    end
    
    line = fgetl(fid);
end
fclose(fid);

cliInfo = struct2table(cliInfo);

% remove patients without survival time
indDel = isnan(cliInfo.time);
cliInfo(indDel, :) = [];

% remove patients with undetermined values
indGrade = strcmp(cliInfo.grade, '[Not Available]') | strcmp(cliInfo.grade, 'GX');
indStage = strcmp(cliInfo.stage, '[Discrepancy]');
indDel = indGrade | indStage;
cliInfo(indDel, :) = [];

% remove patients with time less than 30 days
indDel = cliInfo.time<30;
cliInfo(indDel, :) = [];

save cliInfo cliInfo





