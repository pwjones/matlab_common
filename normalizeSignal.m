function normSig = normalizeSignal(sig)
% function normSig = normalizeSignal(sig)
% Super simple utility, normalizes a vector to a range of 0-1

normSig = sig - nanmin(sig);
normSig = normSig./nanmax(normSig);
