function MiniList = Mickey(Traces)
% Damon Tomlin 12/1/05
%
% Left Mouse Button = add nearby mini
% Right Mouse Button = remove nearby mini
% Left arrow = back one trace
% Right arrow = forward one trace
% Up arrow = forward 100 traces
% Down arrow = backward 100 traces
% q = quit program
% 
% Output is 2D matrix, where is row is:
% [trace# MiniLocation MiniAmplitude]

h = waitbar(0,'Creating Smoothed Traces');
WindowWidth = 20;
S = zeros(size(Traces,2),size(Traces,1));
for i = WindowWidth+1:1:size(Traces,1)-WindowWidth
    S(:,i) = mean(Traces(i-WindowWidth:i+WindowWidth,:),1);
    if (round(i/100)*100 == i)
        waitbar(i/size(Traces,1),h);
    end
end
close(h);    

[StartPoint, EndPoint] = DetermineBoundaries(Traces);
Traces(EndPoint:size(Traces,1),:) = [];
Traces(1:StartPoint,:) = [];
S(:,EndPoint:size(S,2)) = [];
S(:,1:StartPoint) = [];
close(gcf);

h = figure;
MiniList = [];
j = 1;
while j <= size(Traces,2)
    Trace = Traces(:,j);
	Smoothed = S(j,:);
    
    if isempty(MiniList)
        Minis = FindMinis(Smoothed,StartPoint);
    elseif isempty(find(MiniList(:,1) == j))
        Minis = FindMinis(Smoothed,StartPoint);
    else
        Minis = MiniList(find(MiniList(:,1) == j),2:3);
        Minis(find(Minis(:,1) == 0),:) = [];
        MiniList(find(MiniList(:,1) == j),:) = [];
    end
    
    PlotTrace(Trace,Smoothed,Minis,j,size(Traces,2));
    PressType = 0;
    while PressType == 0
        PressType = waitforbuttonpress;
        if PressType == 0
            CP = get(gca,'CurrentPoint');
            x = round(CP(1,1));
            if (x >= 1) & (x <= length(Trace))
                if strcmp(get(gcf,'SelectionType'),'normal')
                    [amp ind] = min(Smoothed(max([x-100 1]):min([x+100 length(Smoothed)])));
                    Minis = [Minis; x+ind-101 amp];
                elseif ~isempty(Minis)
                    Dist = abs(Minis(:,1) - x);
                    Minis(min(find(Dist == min(Dist))),:) = [];
                end
            end
            clf;
            PlotTrace(Trace,Smoothed,Minis,j,size(Traces,2));
        else
            if isempty(Minis)
                MiniList = [MiniList; j 0 0];
            else
                for i = 1:1:size(Minis,1)
                    MiniList = [MiniList; j Minis(i,1) Minis(i,2)];
                end
            end
            MiniList = sortrows(MiniList,[1 2 3]);
            disp(['Trace: ' num2str(j) '   ' num2str(size(Minis,1)) ' Minis']);
            key = double(get(gcf,'CurrentCharacter'));
            switch key
                case 113    %'Q' character quits program
                    j = Inf; 
                case 28     %left arrow moves back 1 trace
                    j = max([j-1 1]);
                case 30     %up arrow moves forward 100 traces
                    j = min([j+100 size(Traces,2)]);
                case 45     %'-' on numpad moves back 10 traces
                    j = max([j-10 1]);
                case 43     %'+' on numpad moves forward 10 traces
                    j = min([j+10 size(Traces,2)]);
                case 31     %down arrow moves back 100 traces
                    j = max([j-100 1]);
                otherwise   %any other character moves forward 1 trace
                    j = min([j+1 size(Traces,2)]);
            end
        end
    end
    clf;
end
MiniList(find(MiniList(:,2) == 0),:) = [];
MiniList(:,2) = MiniList(:,2)-1+StartPoint;


%================================================
function Minis = FindMinis(Smoothed,StartPoint)
i = StartPoint;
D = diff(Smoothed);
Minis = [];
MinNum = 25;
while i<length(D)-MinNum
    if mean(D(i:i+MinNum) < -.25) & sum(D(i:i+MinNum) < -.15) > 20
        [amp, ind] = min(Smoothed(i:min([i+200 length(D)])));
        Minis = [Minis; i+ind amp];
        i = i+ind;
    end
    i = i+1;
end


%================================================
function PlotTrace(Trace,Smoothed,Minis,j,TraceMax)
subplot(2,1,1);
plot(Trace);
hold on;
if ~isempty(Minis)
    plot(Minis(:,1),Minis(:,2),'k.','MarkerSize',25);
end
axis([0 length(Trace) -50 20]);
title(['Trace: ' num2str(j) '/' num2str(TraceMax)],'FontSize',12,'FontWeight','Bold');
subplot(2,1,2);
hold on;
plot(Smoothed);
if ~isempty(Minis)
    plot(Minis(:,1),Minis(:,2),'k.','MarkerSize',25);
end
axis([0 length(Trace) -50 20]);
hold off;


%================================================
function [StartPoint, EndPoint] = DetermineBoundaries(Traces)
SubTraces = [];
for i = 1:100:size(Traces,2)
    SubTraces = [SubTraces; Traces(:,i)'-mean(Traces(:,i))-i];
end
plot(SubTraces');
axis([0 size(Traces,1) -800 50]);
title('Pick Starting Point');
[x,y] = ginput(1);
StartPoint = round(x);
YLim = get(gca,'YLim');
line([StartPoint StartPoint],[YLim(1) YLim(2)],'LineWidth',3);
title('Pick Ending Point');
[x,y] = ginput(1);
EndPoint = round(x);
line([EndPoint EndPoint],[YLim(1) YLim(2)],'LineWidth',3);
drawnow;
clf;


