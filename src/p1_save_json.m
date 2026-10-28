function p1_save_json(name, S)
%P1_SAVE_JSON  Write the numbers of one step to results/<name>.json.
[~, res] = p1_root();
try
    txt = jsonencode(S, 'PrettyPrint', true);
catch
    txt = jsonencode(S);
end
fid = fopen(fullfile(res, [name '.json']), 'w');
fprintf(fid, '%s', txt);
fclose(fid);
fprintf('results written: %s\n', fullfile(res, [name '.json']));
end
