function folders = SetupFolders()

    %Specify required folders
    folders.data = 'Data/';
    folders.sol = 'Data/Solutions/';
    folders.char = 'Data/Characteristics/';
    folders.break = 'Data/Breakage/';
    folders.fig = 'Data/Figures/';


    % Check if there is a data folder
    if ~exist(folders.data, 'dir')
       mkdir(folders.data)
    end

    % Check if there is a solutions folder
    if ~exist(folders.sol, 'dir')
       mkdir(folders.sol)
    end

    % Check if there is a characteristics folder
    if ~exist(folders.char, 'dir')
       mkdir(folders.char)
    end

    % Check if there is a breakage folder
    if ~exist(folders.break, 'dir')
       mkdir(folders.break)
    end

    % Check if there is a figures folder
    if ~exist(folders.fig, 'dir')
       mkdir(folders.fig)
    end




end