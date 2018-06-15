figure;
for ii=100:200
    %polarplot(sniffData.mouseHeadingFromOrtho(ii,:), sniffData.sniffPos(ii,:));
    polarplot(absHeading(ii,:), abs(sniffData.sniffPos(ii,:)));
    hold on;
end