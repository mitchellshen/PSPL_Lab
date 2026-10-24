function pspl_setup()
%PSPL_SETUP  Add the PSPL Lab tool folders to the MATLAB path.
%   Run once per MATLAB session, from any folder:
%
%       run('path/to/PSPL_Lab/pspl_setup.m')
%
%   or type PSPL_SETUP when the repo folder is the current folder.
%   To keep the folders on the path for future sessions, run SAVEPATH
%   afterwards or put the RUN line above in your startup.m.

root = fileparts(mfilename('fullpath'));

% Every folder in the repo, skipping hidden ones such as .git
folders = strsplit(genpath(root), pathsep);
folders = folders(~cellfun(@isempty, folders));
inside  = erase(folders, root);
folders = folders(~contains(inside, [filesep '.']));

addpath(folders{:});
fprintf('PSPL Lab: added %d folders to the MATLAB path.\n', numel(folders));
end
