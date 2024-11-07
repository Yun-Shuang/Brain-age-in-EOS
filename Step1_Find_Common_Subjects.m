
%% Step1: find MSN,FC,SC matrices 
clear;
clc;
%% find common subjects ID for MSN,FC,SC matrices
listpath = '/media/shuang/data/BrainAge/Brain_age_scripts/lists/';
% define lists for patients
file_names = {'FC_list_EOS86.txt', 'MSN_list_EOS95.txt', 'SC_list_EOS90.txt'};
% file_names = {'FC_list_EOS86.txt', 'MSN_list_EOS95.txt', 'PANSS_list_EOS71.txt','SC_list_EOS90.txt'};    
% create a null struct to store all subjects
file_data = cell(1, numel(file_names));
for i = 1:numel(file_names)
    fileName = fullfile(listpath, file_names{i});
    file_data{i} = importdata(fileName);
end
file1_cells = cellstr(file_data{1,1});
file2_cells = cellstr(file_data{1,2});
file3_cells = cellstr(file_data{1,3});
% find common subjects: 80 for patients
common_elements_EOS = intersect(intersect(file1_cells, file2_cells), file3_cells);
writecell(common_elements_EOS,'/media/shuang/data/BrainAge/Glasser360/lists/commonList_EOS.txt')
% file4_cells = cellstr(file_data{1,4});
% common_elements_EOS = intersect(intersect(intersect(file1_cells, file2_cells), file3_cells),file4_cells);
% writecell(common_elements_EOS,'/media/shuang/data/BrainAge/Glasser360/lists/commonList_EOS_PANSS.txt')
% define lists for controls
file_names = {'FC_list_HC91.txt', 'MSN_list_HC98.txt', 'SC_list_HC98.txt'};
file_data = cell(1, numel(file_names));
for i = 1:numel(file_names)
    fileName = fullfile(listpath, file_names{i});
    file_data{i} = importdata(fileName);
end
file1_cells = cellstr(file_data{1,1});
file2_cells = cellstr(file_data{1,2});
file3_cells = cellstr(file_data{1,3});
% find common subjects: 91 for controls
common_elements_HC = intersect(intersect(file1_cells, file2_cells), file3_cells);
writecell(common_elements_HC,'/media/shuang/data/BrainAge/Brain_age_scripts/lists/commonList_HC.txt')
%% estracted common connectivity features for common subjects
% load original connectivity matrices
datapath = '/media/shuang/data/BrainAge/Brain_age_scripts/data/';
load(fullfile(datapath, 'FC_EOS.mat'))
load(fullfile(datapath, 'FC_HC.mat'))
load(fullfile(datapath, 'SC_EOS.mat'))
load(fullfile(datapath, 'SC_HC.mat'))
load(fullfile(datapath, 'MSN_EOS.mat'))
load(fullfile(datapath, 'MSN_HC.mat'))

listpath = '/media/shuang/data/BrainAge/Brain_age_scripts/lists/';
% define lists
file_names = {
    'commonList_EOS.txt',
    'commonList_HC.txt', 
    'FC_list_EOS86.txt',
    'FC_list_HC91.txt',
    'MSN_list_EOS95.txt',
    'MSN_list_HC98.txt',
    'SC_list_EOS90.txt',
    'SC_list_HC98.txt',
    };
% null struct to store subject number
file_data = cell(1, numel(file_names));
for i = 1:numel(file_names)
    file_data{i} = importdata(fullfile(listpath,file_names{i}));
end
commonList_EOS80_cells = (file_data{1,1});
commonList_HC91_cells = cellstr(file_data{1,2});
FC_list_EOS86_cells = cellstr(file_data{1,3});
FC_list_HC91_cells = cellstr(file_data{1,4});
MSN_list_EOS95_cells = cellstr(file_data{1,5});
MSN_list_HC98_cells = cellstr(file_data{1,6});
SC_list_EOS90_cells = cellstr(file_data{1,7});
SC_list_HC98_cells = cellstr(file_data{1,8});
commonList_EOS80_strs = strings(1, length(commonList_EOS80_cells));
for i = 1:length(commonList_EOS80_cells)
    commonList_EOS80_strs(i) = string(commonList_EOS80_cells{i});
end
commonList_HC91_strs = strings(1, length(commonList_HC91_cells));
for i = 1:length(commonList_HC91_cells)
    commonList_HC91_strs(i) = string(commonList_HC91_cells{i});
end
FC_list_EOS86_strs = strings(1, length(FC_list_EOS86_cells));
for i = 1:length(FC_list_EOS86_cells)
    FC_list_EOS86_strs(i) = string(FC_list_EOS86_cells{i});
end
FC_list_HC91_strs = strings(1, length(FC_list_HC91_cells));
for i = 1:length(FC_list_HC91_cells)
    FC_list_HC91_strs(i) = string(FC_list_HC91_cells{i});
end
MSN_list_EOS95_strs = strings(1, length(MSN_list_EOS95_cells));
for i = 1:length(MSN_list_EOS95_cells)
    MSN_list_EOS95_strs(i) = string(MSN_list_EOS95_cells{i});
end
MSN_list_HC98_strs = strings(1, length(MSN_list_HC98_cells));
for i = 1:length(MSN_list_HC98_cells)
    MSN_list_HC98_strs(i) = string(MSN_list_HC98_cells{i});
end
SC_list_EOS93_strs = strings(1, length(SC_list_EOS90_cells));
for i = 1:length(SC_list_EOS90_cells)
    SC_list_EOS93_strs(i) = string(SC_list_EOS90_cells{i});
end
SC_list_HC99_strs = strings(1, length(SC_list_HC98_cells));
for i = 1:length(SC_list_HC98_cells)
    SC_list_HC99_strs(i) = string(SC_list_HC98_cells{i});
end
% store multimodal connectivity features for common subjects
common_FC_HC = zeros(length(commonList_HC91_cells),360,360);
common_SC_HC = zeros(length(commonList_HC91_cells),360,360);
common_MSN_HC = zeros(length(commonList_HC91_cells),360,360);
for i=1:length(commonList_HC91_strs) % HC
    sunj_name = commonList_HC91_strs(i);
    fprintf('Processing file: %s\n', sunj_name);
    idx_in_FC = find(FC_list_HC91_strs == sunj_name);
    idx_in_SC = find(SC_list_HC99_strs == sunj_name);
    idx_in_MSN = find(MSN_list_HC98_strs == sunj_name);

    tmpFC = FCconnectome_HC(:,:,idx_in_FC);
    tmpFC = tmpFC - diag(diag(tmpFC));
    tmpFC = atanh(tmpFC); % Figher's z transformed
    tmpSC = SCconnectome_HC(:,:,idx_in_SC);
    tmpMSN = MSNconnectome_HC(:,:,idx_in_MSN);
    tmpMSN = tmpMSN - diag(diag(tmpMSN));
    tmpMSN = atanh(tmpMSN); % Figher's z transformed

    common_FC_HC(i,:,:) = tmpFC;
    common_SC_HC(i,:,:) = tmpSC;
    common_MSN_HC(i,:,:) = tmpMSN;
end

common_FC_EOS = zeros(length(commonList_EOS80_cells),360,360);
common_SC_EOS = zeros(length(commonList_EOS80_cells),360,360);
common_MSN_EOS = zeros(length(commonList_EOS80_cells),360,360);
for i=1:length(commonList_EOS80_strs) % EOS
    sunj_name = commonList_EOS80_strs(i);
    fprintf('Processing file: %s\n', sunj_name);
    idx_in_FC = find(FC_list_EOS86_strs == sunj_name);
    idx_in_SC = find(SC_list_EOS93_strs == sunj_name);
    idx_in_MSN = find(MSN_list_EOS95_strs == sunj_name);
    
    tmpFC = FCconnectome_EOS(:,:,idx_in_FC);
    tmpFC = tmpFC - diag(diag(tmpFC));
    tmpFC = atanh(tmpFC);
    tmpSC = SCconnectome_EOS(:,:,idx_in_SC);
    tmpMSN = MSNconnectome_EOS(:,:,idx_in_MSN);
    tmpMSN = tmpMSN - diag(diag(tmpMSN));
    tmpMSN = atanh(tmpMSN); % Figher's z transformed

    common_FC_EOS(i,:,:) = tmpFC;
    common_SC_EOS(i,:,:) = tmpSC;
    common_MSN_EOS(i,:,:) = tmpMSN;
end
