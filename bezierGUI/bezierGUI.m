function bezierGUI
	% Bezier GUI
	% To launch, simply type 'play2' in your command window.
	% Alternatively, you can choose run from this editor window, or press 'F5'.
	%
	% To use, choose Add/Delete/Edit from the radio buttons on the right
	% and proceed to edit the bezier curve plot from the GUI. This is done
	% by clicking on the axes. Choose 'Edit' allows to move the points
	% while add adds new points at the points clicked and delete deletes
	% the closest points to your mouse click. You can also change the x and
	% y limits to obtain a better view of the plot.
	%
	% Created on MATLAB 7.10.0 (R2010a) but should work on earlier
	% versions. Makes use of the file Bezier Curve Plotter by Sagar Aiya
	% on the MATLAB File Exchange.
	%
	% (c) Husam Aldahiyat, numandina@gmail.com
	
	
	% Create graphical elements on the GUI
	figure('un','n','pos',[.1 .1 .8 .8],'numbert','off','menub','non',...
		'name','Bezier')
	
	axes('un','n','pos',[.05 .05 .8 .9],'buttondownfcn',@go)
	
	e1 = uicontrol('style','ed','un','n','pos',[.88,.1,.08,.05],...
		'backgroundc',[1,1,1],'string','0 5','fonts',10,'fontn','courier',...
		'callback',@che);
	
	e2 = uicontrol('style','ed','un','n','pos',[.88,.05,.08,.05],...
		'backgroundc',[1,1,1],'string','0 10','fonts',10,'fontn','courier',...
		'callback',@che);
	
	uicontrol('sty','text','un','n','pos',[.86,.045,.02,.05],...
		'string','y','fonts',12,'backgroundc',get(gcf,'color'))
	
	uicontrol('sty','text','un','n','pos',[.86,.095,.02,.05],...
		'string','x','fonts',12,'backgroundc',get(gcf,'color'))
	
	uicontrol('style','text','un','n','pos',[.86,.15,.1,.025],...
		'backgroundc',[1,1,1],'string','Axes Limits','fonts',10,...
		'fontn','courier','backgroundc',get(gcf,'color'));
	
	uicontrol('style','push','un','n','pos',[.86,.375,.1,.05],...
		'backgroundc',[1,1,1],'string','Reset','fonts',10,'fontn','courier',...
		'callback',@gg);
	
	% reset button simply closes and reopens the GUI
	function gg(varargin)
		close(gcf)
		bezierGUI
	end
	
	hg = uibuttongroup('un','n','pos',[.86,.2,.1,.15],...
		'fontn','courier','fonts',10);
	
	h1 = uicontrol('sty','radio','parent',hg,'un','n','pos',[0,1/3,1,1/3],...
		'string','Edit');
	
	h2 = uicontrol('sty','radio','parent',hg,'un','n','pos',[0,2/3,1,1/3],...
		'string','Add');
	
	uicontrol('sty','radio','parent',hg,'un','n','pos',[0,0,1,1/3],...
		'string','Delete');
	
	set(findobj('style','radiobutton'),'fontsize',12)
	
	% changing the x or the y limits
	function che(varargin)
		
		xlim(str2num(get(e1,'string'))); %#ok<*ST2NM>
		ylim(str2num(get(e2,'string')));
		
	end
	
	% original points
	P = [0 0; 1 5; 3 2; 5 10];
	
	% create an initial plot
	bemain(P)
	
	% pressing on the axes
	function go(varargin)
		
		% get points on axes
		xx1 = get(findobj('type','line','marker','o'),'xdata');
		yy1 = get(findobj('type','line','marker','o'),'ydata');
		
		p = double(vpa(get(gca,'currentpoint'),12));
		
		% based on radio button chosen, do something different
		
		% Edit
		if get(hg,'selectedobject') == h1
			
			% get point closest to click
			sst = (abs(p(1,1) - xx1) + abs(p(1,2) - yy1));
			
			% change its x and y to match click point
			P(sst == min(sst),:) = p(1,1:2);
			
		% Add
		elseif get(hg,'selectedobject') == h2
			
			% add point clicked to current points
			P = [P;p(1,1:2)];
			
		% Delete
		else
			
			% get point closest to click
			sst = (abs(p(1,1) - xx1) + abs(p(1,2) - yy1));
			
			% delete it
			P(sst == min(sst),:) = [];
			
		end
		
		% draw the new Bezier curve with the new points
		bemain(P)
		
	end
	
	% taken (comments included) from Bezier Curve Plotter by Sagar Aiya.
	function bemain(P)
		cla
		che
		Px = P(:,1);
		Py = P(:,2);
		
		% Number of points
		r = 100;
		
		% Calculate the Bezier curve
		BezierCurve = bezier(Px,Py,r);
		
		Bx = BezierCurve(:,1);
		By = BezierCurve(:,2);
		
		% Create figure
		
		% Create axes
		grid on
		hold on
		
		% Create plot
		plot(Bx,By,'LineWidth',1) %Bezier curve
		plot(Px,Py,'ro','LineWidth',1); % Control points
		plot(Px,Py,'k-','LineWidth',1); % Lines between control points
		
	end
	
	% taken (comments included) from Bezier Curve Plotter by Sagar Aiya.
	function BezierCurve = bezier(Px,Py,r)
		% This function creates an Bezier Curve with r points and use' the
		% controlpoints Px,Py
		
		% Calculate the degree of the curve
		n = length(Px)-1;
		
		% Calculate the amount of change in the parameter u
		du = 1/r;
		
		% Calculate the binomial koefficienter
		c = zeros(n+1,1);
		bx = zeros(r+1,n+1);
		by = zeros(r+1,n+1);
		
		for i=0:n
			c(i+1,1) = (factorial(n))/(factorial(i)*(factorial(n-i)));
			
			% Calculate the Bernstein polynomium
			for j=0:r
				bx(j+1,i+1) = (c(i+1,1)*(j*du)^(i)*(1-(j*du)).^(n-i))*Px(i+1);
				by(j+1,i+1) = (c(i+1,1)*(j*du)^(i)*(1-(j*du)).^(n-i))*Py(i+1);
			end
		end
		
		Bx = zeros(r+1,1);
		By = zeros(r+1,1);
		% Calculate the Bezier Curve by sum up the Bernstein polynomials
		for k=0:r
			Bx(k+1,1) = sum(bx(k+1,1:n+1));
			By(k+1,1) = sum(by(k+1,1:n+1));
		end
		
		% Output
		BezierCurve = [Bx, By];
	end
	
end