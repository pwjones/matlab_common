function new_areas = mergeAreas(input)
% function new_areas = mergeAreas(input)
% 
% This just creates a new structure with the areas merged. Works for binary areas returned from 'regionprops' with the
% fields 'Area', 'Perimeter', 'PixelList', 'PixelIdxList' so far.


fn = fields(input);
new_areas = input(1);
for jj=1:length(fn)
    if (strcmp(fn(jj), 'Area') || strcmp(fn(jj), 'Perimeter'))
        new_areas.(fn{jj}) = sum([input.(fn{jj})]);
    elseif (strcmp(fn(jj), 'PixelIdxList') || strcmp(fn(jj), 'PixelList'))
        temp = [];
        for ii = 1:length(input)
            temp = cat(1, temp, input(ii).(fn{jj}));
        end
        new_areas.(fn{jj}) = temp;
    end
end
