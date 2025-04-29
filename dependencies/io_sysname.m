% io_sysname.m
% Zoelch Archibald, 2025.
%
% USAGE:
% [metname]=io_sysname({strings});
% 
% DESCRIPTION:
% Returns name of the metabolite stored in sys.name. As sys can be a
% structure array {sys.name} can contain several names. In this case, this functions assumes that the
% name is of the form met_(some description of this part of the met). So
% this functions reduces the name to met.
% 
% INPUTS:
% {strings}   = cell array of strings
%
% OUTPUTS:
% metname = name of the metabolite.


function metname = io_sysname(strings)
% 
    % wenn nur 1 string dann einfach der name der gegeben ist
    % wenn es mehrere hat, dann str_split
    % das heisst man muss beim naming darauf achten-> dort hinschreiben
    % Absichern, dass es eine Zelle von Strings/Chars ist
    if ~iscell(strings)
        error('Eingabe muss eine Zell-Array von Strings oder Char-Vektoren sein.');
    end

    % Leere oder nur ein Element?
    if isempty(strings)
        metname = '';
        return;
    elseif numel(strings) == 1
        metname = strings{1};
        return;
    end

    % Sonst
    % Starte mit dem ersten String als Referenz
    first_name_split=strsplit(strings{1},'_');
    metname =first_name_split{1};
end
