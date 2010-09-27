function cell_csv_write(filename, C)
fid=fopen(filename, 'w');
s=size(C)
for x=1:s(1)
    for y=1:s(2)
        if ischar(C{x,y})
            if strcmp(C{x,y},'')
                fprintf(fid, '; ');
            else
                fprintf(fid, '"%s"; ', C{x,y});
            end
        elseif isinteger(C{x,y})
            fprintf(fid, '%15d; ', C{x,y});
        elseif isfloat(C{x,y})
            fprintf(fid, '%15.5f; ', C{x,y});
        else
            fprintf(fid, '; ');
        end
    end
    fprintf(fid, '\n');
end
fclose(fid);