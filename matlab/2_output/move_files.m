% Comment
subDirs = {'Version v1 (72.2)' 'Version v2 (72.2)' 'Version v3 (72.2)'};

disp(mfilename('fulllpath'))

files = dir(fullfile(fileparts(fileparts(pwd)), '**/*.mat'));
files = files(~[files.isdir]);

for i = 1:length(files)
    load(files(i).name), 'PROBLEM_CONSTANTS)';
    fname = func2str(PROBLEM_CONSTANTS.function_to_minimize);
    script_name = files(i).name;
    
    if contains(fname, 'v1')
        movefile(fullfile(pwd, script_name), fullfile(pwd, subDirs{1}, script_name));
    end
     if contains(fname, 'v2')
        movefile(fullfile(pwd, script_name), fullfile(pwd, subDirs{2}, script_name));
     end
     if contains(fname, 'v3')
        movefile(fullfile(pwd, script_name), fullfile(pwd, subDirs{3}, script_name));
    end
end