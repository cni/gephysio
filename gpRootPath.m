function rootPath = gpRootPath()
% Return the path to the gephysio directory
%
% Syntax:
%    rootPath = gpRootPath;
%
% Description:
%    Return base bath of the gephysio directory.
%
% Inputs:
%    None.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% See Also:

%%  chdir(fullfile(gpRootPath,'local'))

%%
rootPath = fileparts(which('gpRootPath'));


end