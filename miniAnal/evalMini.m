function f = evalMini(t, timeVect, miniVect)
% function f = evalMini(x, timeVect, miniVect)
%
% Function evaluates a mini of given miniVect,at a certain time t, 
% where the index of time t is determined by timeVect, a time vector the 
% same size as miniVect
f=NaN*ones(size(t));
if(size(timeVect) ~= size(miniVect))
    disp('evalMini: Vect sizes not equal');
    return;
end
for i =1:length(t)
    iv = find(timeVect >= t(i)); 
    if (isempty(iv))
        %disp('evalMini: t is not contained in the given timeVect.');
        f(i) = NaN;
    else
        i_match = iv(1);
        f(i) = miniVect(i_match);
    end
end