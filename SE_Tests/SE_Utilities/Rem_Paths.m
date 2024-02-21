function Rem_Paths(additional_paths)

for k=1:length(additional_paths)
    rmpath(additional_paths{k});   
end