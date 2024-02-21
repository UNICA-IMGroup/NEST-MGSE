function additional_paths = Add_Paths()

file_folder='.\';
additional_paths = {file_folder, ...
    '..\SE_Library\', '..\SE_Library\Core\', '..\SE_Library\Networks\', '..\SE_Library\Measurements\', '.\SE_Utilities\','.\TestsConfig\', ...   
    }; 

for k=1:length(additional_paths)
    addpath(additional_paths{k});   
end