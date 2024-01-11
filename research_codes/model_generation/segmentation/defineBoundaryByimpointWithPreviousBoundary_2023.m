% clear all; close all; clc;
% updated by Arash

function [xycoor_, hpts]=defineBoundaryByimpointWithPreviousBoundary_2023(imData, bc_data, closeB, str_title, lv_p)
global finish xycoor xcoor ycoor
finish=0;
refresh=0;

fig = figure;
imshow(imData,[]); hold on;
title(str_title);

%%then we can plot bc_data 
if ~isempty(bc_data)
    hold on;
    plot(bc_data(1,:), bc_data(2,:), 'r-', 'LineWidth', 2);
end

if ~isempty(lv_p)
    hold on;
    plot(lv_p(1,:), lv_p(2,:), 'b--', 'LineWidth', 1);
end

pos_fig = get(gcf,'Position');
set(0,'units','pixels');Pix_SS =get(0,'screensize');


%imSize = size(im);
% pos_tb1 = [20, pos_fig(2)+pos_fig(4)/2, 60, 20];
% pos_tb1 = [pos_fig(1)+pos_fig(3)/2-105, pos_fig(2)+pos_fig(4)+30, 90, 30];
pos_tb1 = [Pix_SS(3)/2-105, Pix_SS(4)-150, 90, 30];


tb1 = uicontrol(fig, 'Style', 'togglebutton', 'String', 'Stop', ...
                 'Position', pos_tb1, 'BackgroundColor','g','FontSize',12);
pos_tb1 = tb1.Position;
pos_tb2 = [pos_tb1(1)+30+pos_tb1(3), pos_tb1(2), pos_tb1(3),  pos_tb1(4)];

tb2 = uicontrol(fig, 'Style', 'pushbutton', 'String', 're-plot', ...
                'Position', pos_tb2, 'BackgroundColor','g',...
                'Callback', @replot,'FontSize',12);

            
            
            
set(gcf,'units','normalized','outerposition',[0 0 1 1]);            
set(gcf, 'Windowkeypressfcn',@keyboard_shortcuts);
set(gcf, 'WindowButtonUpFcn',@mouse_shortcuts)
pIndex = 0;
while tb1.Value==0
  pIndex = pIndex + 1;  
  hpts(pIndex).h = impoint(gca);
end

xycoor = zeros([2, pIndex]);

figData.ptlist = hpts;
figData.hfig = fig;
figData.pIndex = pIndex;
%figData.hcurve = hfig;

%%press anykey twice to return to the main function
%       w = waitforbuttonpress;
      while true
      if finish==1
          break;
      end
%       if tb1.Value==1
%           replot; drawnow;
%       end
waitforbuttonpress;
      end
      
%       %%update the xycoor
      xycoor(1,:) = xcoor';
      xycoor(2,:) = ycoor';
      xycoor_=xycoor;
%       hpts = hptsT;
hpts =[];
      close(fig);

    % Create Nested Callback Functions (Programmatic Apps)
    function replot(hObject,eventdata)
refresh=1;
        %%delete the old curve
        if isfield(figData, 'hcurve')
            delete(figData.hcurve);
        end
         
        %plot the new curve
        hptsT = figData.ptlist ;
        pIndexT = figData.pIndex;
        
        xcoor = [];
        ycoor = [];
        
        for i = 1 : pIndexT
            pos = getPosition(hptsT(i).h);
            xcoor(i,1) = pos(1);
            ycoor(i,1) = pos(2);
        end
        
        %%reinterpolate with spline to have a better curve
        if closeB 
            tt = 1:0.1:pIndexT+1;
            t = 1 : pIndexT+1;
            xx = spline(t, [xcoor' xcoor(1)], tt);
            yy = spline(t, [ycoor', ycoor(1)], tt);
        else
            tt = 1:0.1:pIndexT;
            t = 1 : pIndexT;
            xx = spline(t, xcoor', tt);
            yy = spline(t, ycoor', tt);
            
        end

        figure(figData.hfig); hold on;
        %hcurve = plot([xcoor' xcoor(1)], [ycoor', ycoor(1)], 'g-', ...
        %    'LineWidth',2);
        
        hcurve = plot(xx, yy, 'g-', ...
            'LineWidth',.5);

        % for i = 1 : pIndex
        %     uistack(hpts(i).h, 'top');
        % end


         uistack(hcurve,'bottom');
         uistack(hcurve,'up',1);
         figData.hcurve = hcurve;
         figData.ptlist = hptsT;
        % uistack(fig, 'bottom');
        
        
        %%update the xycoor
        xycoor(1,:) = xcoor';
        xycoor(2,:) = ycoor';
       

        
    end
        function keyboard_shortcuts(source,eventdata)
%             disp(eventdata)
keyPressed = eventdata.Key;
      if strcmpi(keyPressed,'r')
          replot;
      end
            if strcmpi(keyPressed,'s')
          tb1.Value=1;
            end
          if strcmpi(keyPressed,'space')
finish=1;      
      end
        end
 function mouse_shortcuts(source,eventdata)
            
      if tb1.Value==1 && refresh
          
          replot;
      end

end

end






%  delete(hcurve)


