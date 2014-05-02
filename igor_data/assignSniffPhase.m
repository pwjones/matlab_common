function phase_vect = assignSniffPhase(inhaleVect)
% function phase = assignSniffPhase(inhaleVect)
%
% Simple function - given a vector of zeros and 1s to denote sniff
% times (inhalation peak is zero).

phase_vect = NaN*zeros(size(inhaleVect));
inhale = find(inhaleVect);
for ii = 1:(length(inhale)-1) 
    in1 = inhale(ii);
    in2 = inhale(ii+1);
    tp = linspace(0,2*pi, in2-in1+1);
    phase_vect(in1:(in2-1)) = tp(1:end-1);
end
