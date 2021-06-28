function call_tetgen(filename, refinement)
%CALL_TETGEN Call Tetgen executable from system with refinement.
%   The tetgen executable is chosen according to the operating system. If an
%   error occurs, you might have to manually change the rights of the Tetgen
%   executable.
%
%   filename: string
%   refinement: [1 x 1]


% Tetgen command for corresponding operating system

tmpstr=which('call_tetgen');
tmpstr=tmpstr(1:end-22);


if ispc
    tetgen_cmd = fullfile(tmpstr,"tetgen\win64\tetgen");
elseif ismac
    tetgen_cmd = fullfile(tmpstr,"tetgen/mac64/tetgen");
elseif isunix
    tetgen_cmd = fullfile(tmpstr,"tetgen/lin64/tetgen");
else
    warning("Using Linux Tetgen command.")
    tetgen_cmd = fullfile(tmpstr,"tetgen/lin64/tetgen");
end

% Options for Tetgen command
if  nargin == nargin(@call_tetgen) && refinement > 0
    % Pass refinement to the 'a' flag of Tetgen. This gives a maximum
    % tetrahedron volume (not length, as in earlier versions)
    % tetgen_options = "-pq0.5AVa" + num2str(setup.pde.refinement);
    tetgen_options = "-pqAVa" + num2str(refinement);
else
    tetgen_options = "-pqAV";
end

% Call Tetgen
cmd = sprintf("%s %s %s", tetgen_cmd, tetgen_options, filename);
disp(cmd)
system(cmd);
