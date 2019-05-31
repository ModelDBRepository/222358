function [dPatch,ps] = draggableRect(varargin)
% Create a draggable rectangle; optionally snap to a grid
%
% NOTE: All positions are specified in NORMALIZED units.
% (Unlike MATLAB's DRAGRECT, which does not work in
% normalized units.)
%
% SYNTAX:
% draggableRect
%    Creates a draggable, unconstrained (no snapTo)
%    rectangle in the current axes. (Drag to define the
%    initial position.)
%
% [dRect,ps] = draggableRect
%    Also returns a handle, dRect, to the rectangle, and an
%    array of handles, ps, to the draggable points. ps are
%    in the order: [BL,BR,TL,TR,BM,TM,RM,LM].
%    (B=Bottom;R=Right;T=Top;L=Left;M=Middle).
%
% [...] = draggableRect(...'snapTo',TF)
%    snapTo: a logical bit (T/F, or 1/0) indicating whether
%       rectangle should snap to a normalized grid, with x-
%       and y- gridpoints specified by [0:0.1:1]. Default =
%       false.
%
% [...] = draggableRect(...,'xs',xs,'ys',ys)
%       If snapTo is True, you may also specify additional
%       arguments xs and ys, vectors of x-positions and
%       y-positions, respectively, to snap to. All positions
%       should be specified in normalized units. (Note that
%       you may also specify only a single vector, which
%       will be used for both xs and ys.)
%
% [...] = draggableRect(...,'parent',parent)
%       Optionally, you may specify the handle to the axis
%       in which you want to create the draggable rectangle.
%       This argument can go anywhere in the inputs.
%
% [...] = draggableRect(...,'InitialPosition',[x y w h])
%       Optionally, you may specify the handle to the axis
%       in which you want to create the draggable rectangle.
%       This argument can go anywhere in the inputs.
%
% [...] = draggableRect(...,'newPositionCallback',fcnHndl)
%       Optionally, you may specify a function handle, fcnHndl, which will
%       be executed when the rectangle is dragged.
%
% EXAMPLES:
%
%%% Ex.1
%%% Create an unconstrained draggable rectangle in the specified axes.
%
% gridAx = axes('units','normalized','position',[0.05 0.05 0.9 0.9],...
%                'xlimmode','manual','ylimmode','manual');
% draggableRect(false,gridAx);
%
%%% Ex.2
%%% Create a constrained (snapTo-on) draggable rectangle
%
% gridAx = axes('units','normalized','position',[0.05 0.05 0.9 0.9]);
% % Create grid array
% [gridXs,gridYs] = meshgrid(0:20,20:-1:0);
% % Normalize by dividing by maximum:
% gridXs = gridXs/20;
% gridYs = gridYs/20;
% % Visualize grid array:
% plot(gridXs,gridYs,'linestyle','none','marker','o',...
%    'markersize',4,'hittest','off');
% [r,p] = draggableRect('snapTo',true,'xs',gridXs(1,:),'ys',gridYs(:,1));
%
%%% Ex.3
%%% Trigger a function callback when dragged 
%
% gridAx = axes('units','normalized','position',[0.05 0.05 0.9 0.9],...
%                'xlimmode','manual','ylimmode','manual');
% cbFcn = @(x) title(sprintf('Lower Left: (%0.1f, %0.1f)',x(1),x(2)))
% draggableRect(false,gridAx,'initialPosition',[0.1 0.1 0.2 0.2],...
%                'newPositionCallback',cbFcn);
%
%
%%% Ex.4
%%% Drag a rectangle top of an image
%%% Note that the draggableRect must be in a [0,1]-limited axes; to use it
%%%    on top of an image (for example), create a separate axes on top of
%%%    the image-containing one:
%
% figure;
% imshow('peppers.png');
% imax = gca;
% axes('parent',gcf,... 
% 'units','normalized',... 
% 'position',get(imax,'position'),... 
% 'xlimmode','manual','ylimmode','manual',... 
% 'xtick',[],'ytick',[],'color','none'); 
% draggableRect('initialPosition',[0.1 0.1 0.2 0.2]);
%
%
%%
% NOTE: If you want to work with NON-NORMALIZED units, use
% MATLAB's |dragrect| instead.
%
% See also: dragrect

% Written by Brett Shoelson, PhD
% brett.shoelson@mathworks.com
% 1/10/2012
% Acknowledgement: Francois Bouffard, for his excellent
%    implementation of <draggable>.
%
% Modifications:
% 5/28/2014: Now allows specification of a newPositionCallbackFunction.
%    (Thanks to Eric Deal for the prompt.) Added two new examples--one to
%    show a callback function, one to show how to use draggableRect on top
%    of an image.
%
% Copyright 2014 The MathWorks, Inc.

[gridAx,snapTo,xs,ys,initialPos,fcnHndl] = parseInputs(varargin);
% Define a context menu; it is not attached to anything
hcmenu = uicontextmenu;
% Define the context menu items and install their callbacks
opts = {'off','on'};
uimenu(hcmenu, 'Label', 'Snap to Grid', 'Callback', @toggleSnap,'checked',opts{snapTo+1});
uimenu(hcmenu, 'Label', 'Get Position', 'Callback', @getRectPos);
uimenu(hcmenu, 'Label', 'Set Position', 'Callback', @setRectPos);
uimenu(hcmenu, 'Label', 'Lock Position', 'Callback', @lockIt);
uimenu(hcmenu, 'Label', 'Lock Aspect Ratio', 'Callback', @lockAspectRatio);
uimenu(hcmenu, 'Label', 'Change Rectangle Display',  'Callback', @inspectRect);
uimenu(hcmenu, 'Label', 'Delete',  'Callback', @deleteIt);

parentFig = ancestor(gridAx,'figure');
cPointer = get(parentFig,'pointer');
if snapTo
    dx = xs(2)-xs(1);
    dy = ys(2)-ys(1);
    %point1 = snap(point1,xs,ys);
else
    %Defaults for minimum rectangle size
    dx = 0.02;
    dy = 0.02;
end
if isempty(initialPos)
set(parentFig,'pointer','crosshair');
waitforbuttonpress;
point1 = get(gridAx,'CurrentPoint');    % button down detected
if snapTo
    point1 = snap(point1,xs,ys);
end
rbbox;                               % return figure units
set(parentFig,'pointer',cPointer);
point2 = get(gridAx,'CurrentPoint'); % button up detected
if snapTo
    point2 = snap(point2,xs,ys);
else
    xs = [];
    ys = [];
end
point1 = point1(1,1:2);              % extract x and y
point2 = point2(1,1:2);
dim = abs(point1-point2);            % and dimensions
p1 = min(point1,point2);             % calculate locations
p = [p1(1) p1(2) dim(1) dim(2)];
p = max(p,[0 0 abs(dx),abs(dy)]);
else
    p = initialPos;
end
r = p2r(p);

% CREATE PATCH...
dPatch = patch(r(1,:),r(2,:),[0.59 1 0.875],...
    'facealpha',0.5,'tag','draggablePatch');
% ...and make it draggable
draggable(dPatch,@dragged)
snapCtrls.snapTo = snapTo;
snapCtrls.xs = xs;
snapCtrls.ys = ys;
snapCtrls.dx = dx;
snapCtrls.dy = dy;
% Store snap-to and active-state information
setappdata(dPatch,'snapCtrls',snapCtrls);
setappdata(dPatch,'isActive',true)
%
edgePoints = [];
centroid = [];
% Calculate the positions of the edgePoints:
positionEdgePoints(p);
% Create edgePoints, store handles in patch's appdata:
drawEdgePoints;
%
set(dPatch,'uicontextmenu',hcmenu)

    function deleteIt(varargin)
        delete(dPatch,ps);
    end

    function dragged(varargin)
        %[BL,BR,TL,TR,BM,TM,RM,LM]
        C1 = getappdata(dPatch,'Centroid');
        snapCtrls = getappdata(dPatch,'snapCtrls');
        snapTo = snapCtrls.snapTo;
        xs = snapCtrls.xs;
        ys = snapCtrls.ys;
        dx = snapCtrls.dx;
        dy = snapCtrls.dy;
        pt = varargin{1};
        handleType = get(pt,'tag');
        if strcmp(handleType,'draggablePatch')
            %oldp = p;
            r = [get(pt,'xdata'),get(pt,'ydata')]';
            p = r2p(r);
            if snapTo
                newxy = snap([p(1),p(2)],xs,ys);
                p(1) = newxy(1);
                p(2) = newxy(2);
                r = p2r(p);
                set(pt,'xdata',r(1,:),'ydata',r(2,:));
            end
            p = min(p,[1-p(3) 1-p(4) Inf Inf]);
            r = p2r(p);
            positionEdgePoints(p);
        else
            newx = get(pt,'xdata');
            newy = get(pt,'ydata');
            if snapTo
                newxy = snap([newx,newy],xs,ys);
                newx = newxy(1);
                newy = newxy(2);
            end
            r = [get(dPatch,'xdata'),get(dPatch,'ydata')]';
            cpos = r2p(r);
            Lx = cpos(1);
            By = cpos(2);
            dim = cpos(3:4);
            Rx = dim(1)+Lx;
            Ty = dim(2)+By;
            switch handleType
                case 'BL'
                    Lx = max(0,min(newx,Rx-abs(dx)));
                    By = max(0,min(newy,Ty-abs(dy)));
                case 'BR'
                    Rx = min(1,max(newx,Lx+abs(dx)));
                    By = max(0,min(newy,Ty-abs(dy)));
                case 'TL'
                    Lx = max(0,min(newx,Rx-abs(dx)));
                    Ty = min(1,max(newy,By+abs(dy)));
                case 'TR'
                    Rx = min(1,max(newx,Lx+abs(dx)));
                    Ty = min(1,max(newy,By+abs(dy)));
                case 'BM'
                    By = max(0,min(newy,Ty-abs(dy)));
                case 'TM'
                    Ty = min(1,max(newy,By+abs(dy)));
                case 'RM'
                    Rx = min(1,max(newx,Lx+abs(dx)));
                case 'LM'
                    Lx = max(0,min(newx,Rx-abs(dx)));
            end
            dim = [Rx-Lx Ty-By];
            p = [Lx By dim(1) dim(2)];
            p = max(p,[0 0 abs(dx),abs(dy)]);
            positionEdgePoints(p);
            r = p2r(p);
        end
        set(dPatch,'xdata',r(1,:),'ydata',r(2,:));
        redrawEdgePoints;
        C2 = getappdata(dPatch,'Centroid');
        setappdata(dPatch,'del',[C2(1)-C1(1),C2(2)-C1(2)]);
		if ~isempty(fcnHndl)
			fcnHndl(p)
		end
    end

    function drawEdgePoints
        % Creates edgePoints and stores them in patch's
        % appdata
        oldhold = ishold;
        hold on
        ps(1) = plot(edgePoints.BL(1),edgePoints.BL(2),'bs','tag','BL');
        draggable(ps(1),@dragged);
        ps(2) = plot(edgePoints.BR(1),edgePoints.BR(2),'bs','tag','BR');
        draggable(ps(2),@dragged);
        ps(3) = plot(edgePoints.TL(1),edgePoints.TL(2),'bs','tag','TL');
        draggable(ps(3),@dragged);
        ps(4) = plot(edgePoints.TR(1),edgePoints.TR(2),'bs','tag','TR');
        draggable(ps(4),@dragged);
        ps(5) = plot(edgePoints.BM(1),edgePoints.BM(2),'bs','tag','BM');
        draggable(ps(5),'constraint','vertical',@dragged);
        ps(6) = plot(edgePoints.TM(1),edgePoints.TM(2),'bs','tag','TM');
        draggable(ps(6),'constraint','vertical',@dragged);
        ps(7) = plot(edgePoints.RM(1),edgePoints.RM(2),'bs','tag','RM');
        draggable(ps(7),'constraint','horizontal',@dragged);
        ps(8) = plot(edgePoints.LM(1),edgePoints.LM(2),'bs','tag','LM');
        draggable(ps(8),'constraint','horizontal',@dragged);
        set(ps,'busyaction','cancel');
        setappdata(dPatch,'edgePoints',ps)
        if ~oldhold
            hold off
        end
    end

    function getRectPos(varargin)
        disp(p)
        disp('Position p written to base workspace')
        assignin('base','p',p)
    end

    function inspectRect(varargin)
        inspect(gco)
    end

    function lockAspectRatio(varargin)
        %[BL,BR,TL,TR,BM,TM,RM,LM]
        menuItem = varargin{1};
        currARLocked =  get(menuItem,'checked');%uicontextmenu
        if strcmp(currARLocked,'off')
            angle = atan(p(4)/p(3));
            set(menuItem,'checked','on')
            set(ps(5:end),'visible','off')
            draggable(ps([1,4]),@dragged,'constraint','diagonal',angle);
            draggable(ps(2:3),@dragged,'constraint','diagonal',-angle);
        else
            set(menuItem,'checked','off')
            set(ps,'visible','on')
            draggable(ps(1:4),@dragged,'constraint','none');
        end
    end

    function lockIt(varargin)
        menuItem = varargin{1};
        currLock =  get(menuItem,'checked');%uicontextmenu
        if strcmp(currLock,'off')
            set(menuItem,'checked','on')
            draggable(ps,'off')
            set(ps,'color',[0.7 0.7 0.7])
            draggable(dPatch,'off')
        else
            set(menuItem,'checked','off')
            draggable(dPatch,@dragged)
            draggable(ps(1:4),@dragged);
            draggable(ps(5:6),'constraint','vertical',@dragged);
            draggable(ps(7:8),'constraint','horizontal',@dragged);
            set(ps,'color','b')
        end
    end

    function positionEdgePoints(p)
        % Calculates new values for edgePoint positions
        dim = [p(3) p(4)];
        Lx = p(1);
        Rx = dim(1)+Lx;
        Mx = 0.5*(Lx+Rx);
        By = p(2);
        Ty = dim(2)+By;
        My = 0.5*(Ty+By);
        %
        edgePoints.BL = [Lx, By];%Bottom Left
        edgePoints.BR = [Rx, By];%Bottom Right
        edgePoints.TL = [Lx, Ty];%Top Left
        edgePoints.TR = [Rx, Ty];%Top Right
        %
        edgePoints.BM = [Mx, By];%Bottom Middle
        edgePoints.RM = [Rx, My];%Right Middle
        edgePoints.TM = [Mx, Ty];%Top Middle
        edgePoints.LM = [Lx, My];%Left Middle
        %
        centroid =      [Mx, My];
        setappdata(dPatch,'Centroid',centroid);
    end

    function redrawEdgePoints
        % Gets calculated edgePoint positions and moves points accordingly
        Lx = edgePoints.BL(1);
        By = edgePoints.BL(2);
        Rx = edgePoints.BR(1);
        Ty = edgePoints.TL(2);
        Mx = edgePoints.BM(1);
        My = edgePoints.LM(2);
        %[BL,BR,TL,TR,BM,TM,RM,LM]
        set(ps(1),'xdata',Lx,'ydata',By);
        set(ps(2),'xdata',Rx,'ydata',By);
        set(ps(3),'xdata',Lx,'ydata',Ty);
        set(ps(4),'xdata',Rx,'ydata',Ty);
        set(ps(5),'xdata',Mx,'ydata',By);
        set(ps(6),'xdata',Mx,'ydata',Ty);
        set(ps(7),'xdata',Rx,'ydata',My);
        set(ps(8),'xdata',Lx,'ydata',My);
    end

    function setRectPos(varargin)
        p = promptForPosition(p);
        r = p2r(p);
        set(dPatch,'xdata',r(1,:),'ydata',r(2,:));
        positionEdgePoints(p);
        redrawEdgePoints
    end

    function toggleSnap(varargin)
        currPatch = gco; %rect
        currSnap =  get(varargin{1},'checked');%uicontextmenu
        snapCtrls = getappdata(currPatch,'snapCtrls');
        if strcmp(currSnap,'off')
            set(varargin{1},'checked','on')
            snapCtrls.snapTo = true;
        else
            set(varargin{1},'checked','off')
            snapCtrls.snapTo = false;
        end
        setappdata(currPatch,'snapCtrls',snapCtrls);
    end
end

function [p,palt,dx,dy] = snap(p,xs,ys)
% Note: I use MIN here because in the calculation of the
% alternate output I set a value to infinity. If that
% happens to be in the 1st or 2nd location, that would cause
% a problem. And I pass out dx and dy (without using them)
% so that the version in the caller workspace is not
% modified.
dx = min(xs(2)-xs(1),xs(4)-xs(3));
dy = min(ys(2)-ys(1),ys(4)-ys(3));
if size(p,1) == 1 %edgePoint drag triggers the resize
    %              (vs instantiation via get(gca,'CurrentPoint'))
    [~,indx] = min(abs(xs-p(1)));
    [~,indy] = min(abs(ys-p(2)));
else
    [~,indx] = min(abs(xs-p(2,1)));
    [~,indy] = min(abs(ys-p(2,2)));
end
p(1,1) = xs(indx);
p(1,2) = ys(indy);
if nargout > 1
    % Second-closest requested
    % NOTE: Because I pass in xs and ys as inputs, modifying
    % them here shouldn't affect them in the global scheme.
    xs(indx) = Inf;
    ys(indy) = Inf;
    palt = snap(p,xs,ys);
end
end

function [ax,snapTo,xs,ys,initialPos,fcnHndl] = parseInputs(varargin)
inputs = varargin{:};
% Defaults
ax = gca;
snapTo = false;
xs = 0:0.25:1;
ys = 0:0.25:1;
initialPos = [];
fcnHndl = [];
for ii = 1:2:numel(inputs)-1 
    switch lower(inputs{ii})
        case {'ax','axes'}
            ax = inputs{ii+1};
        case 'snapto'
            snapTo = inputs{ii+1};
        case 'xs'
            xs = inputs{ii+1};
        case 'ys'
            ys = inputs{ii+1};
        case {'initpos','initialpos','initialposition','initposition'}
            initialPos = inputs{ii+1};
            validPos = @(x) all([all(isnumeric(x)), numel(x)==4]);
            if ~validPos(initialPos)
                error('DRAGGABLERECT: Initial position must be specified in normalized values');
			end
		case 'newpositioncallback'
			fcnHndl = inputs{ii+1};
    end
end
end

function r = p2r(p)
r = [p(1) p(1)+p(3) p(1)+p(3) p(1) p(1);
    p(2) p(2) p(2)+p(4) p(2)+p(4) p(2)];
end

function  p = promptForPosition(varargin)
p = varargin{1};
oldP = p;
prompt={'x','y','w','h'};
name='Define Normalized Position: [x y width height]';
defaultanswer={num2str(p(1)),num2str(p(2)),num2str(p(3)),num2str(p(4))};
tmp=inputdlg(prompt,name,1,defaultanswer);
if ~isempty(tmp)
    p(1) = str2double(tmp{1});
    p(2) = str2double(tmp{2});
    p(3) = str2double(tmp{3});
    p(4) = str2double(tmp{4});
end
if any(isnan(p)) || any(p < 0) || any(p > 1)
    tmp = questdlg('At least one value you specified is not on [0,1]. Do you want to continue?','WARNING!','CANCEL','Continue','CANCEL');
    if strcmp(tmp,'CANCEL')
        p = oldP;
    end
end
end

function p = r2p(r)
p = [r(1,1) r(2,1) r(1,2)-r(1) r(2,3)-r(2,2)];
end