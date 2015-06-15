%%This Script Simulates a double pendulum, using lagrangian calculations
close all;
% clc;

%Settings
save = false;       %save to .avi?
skip=2;           %playback speed.

q = get(simOut, 'q');
q = squeeze(q.signals.values);
signalSize = size(q);
rows = max(signalSize);
signalSize = min(signalSize)/2;
the = q(:,1:signalSize/2);
r =   q(:,signalSize/2+1:signalSize);


%% Computing the xy coordinates of the bobs for my own solution
x = zeros(rows,signalSize/2);
y = zeros(rows,signalSize/2);
xm = zeros(rows,signalSize/2);
ym = zeros(rows,signalSize/2);

for i = 1:signalSize/2
    x(:,i) = r(:,i).*cos(sum(the(:,1:i),2));
    y(:,i) = r(:,i).*sin(sum(the(:,1:i),2));
    xm(:,i) = x(:,i)/2;
    ym(:,i) = y(:,i)/2;
end

for i = 2:signalSize/2
    x(:,i) = sum(x(:,1:i),2);
    y(:,i) = sum(y(:,1:i),2);
end

for i = 2:signalSize/2
    xm(:,i) = xm(:,i) + x(:,i-1);
    ym(:,i) = ym(:,i) + y(:,i-1);
end


% yx = y.signals.values(:,1);
% yy = y.signals.values(:,2);

% u1 = sin(u(:,1)).*l;
% v1 = -cos(u(:,1)).*l;
% u2 = u1 + sin(u(:,1)+u(:,2)).*l;
% v2 = v1 - cos(u(:,1)+u(:,2)).*l;

%% Plotting the response in xy
F(rows) = struct('cdata', [ ], 'colormap', [ ]);
hFig = figure('Name','Anim');
if(save)
    set(hFig,'Visible','off');
    aviobj=VideoWriter('animation');
    open(aviobj);
end

% set(ax,'nextplot','replacechildren');
% set(ax, 'XLimMode', 'manual');
% set(ax, 'YLimMode', 'manual');


xymin = -4;
xymax = 4;
F(max(size(1:skip:size(the,1)))) = struct('cdata',[],'colormap',[]);
if(save)
    for k = 1:skip:max(size(the));
        plot([0 x(k,:)] , [0 y(k,:)],'-o');
        hold on;
        plot(xm(k,:), ym(k,:),'ok');
%         plot(yx(k), yy(k),'or');
        hold off;

        axis([xymin, xymax, xymin, xymax]);
        xlabel('x[m]');
        ylabel('y[m]');

        img = hardcopy(hFig, '-dzbuffer', '-r0');
        writeVideo(aviobj, im2frame(img));
    end
else
    for k = 1:skip:max(size(the));
        plot([0 x(k,:)] , [0 y(k,:)],'-o');
        hold on;
        plot(xm(k,:), ym(k,:),'ok');
%         plot(yx(k), yy(k),'or');
        hold off;

        axis([xymin, xymax, xymin, xymax]);
        xlabel('x[m]');
        ylabel('y[m]');

        axis([xymin, xymax, xymin, xymax]);
        xlabel('x[m]');
        ylabel('y[m]');

        F(k) = getframe(hFig);
    end
end


if(save)
    close(aviobj);
end
    %% Playing the animation
    %movie (F, 1 , 25);