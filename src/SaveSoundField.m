function SaveSoundField(filename, tlmin, tlmax, r, z, tl)

    tl_fid = fopen(filename, 'w');

    fwrite(tl_fid, length(z), 'int32');
    fwrite(tl_fid, length(r), 'int32');
    fwrite(tl_fid, tlmin, 'double');
    fwrite(tl_fid, tlmax, 'double');
    fwrite(tl_fid, z,     'double');
    fwrite(tl_fid, r,     'double');
    fwrite(tl_fid, tl,    'double');
    
    fclose(tl_fid);
end
