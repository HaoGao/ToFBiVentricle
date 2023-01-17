function [endoT, epiT]= manualSegEndoEpiNoSampling(im)

endoT = [];
h=figure(); hold on; imshow(im,[]);hold on;
 but = 1; n = 1; endo=[];
 while but == 1
                [x, y, but] = ginput(1);
                plot(x,y,'r*');
                endo(n,1)=x;
                endo(n,2)=y;
                n = n+1;
 end
 
 but = 1; n = 1; epi=[];
 while but == 1
                [x, y, but] = ginput(1);
                plot(x,y,'b*');
                epi(n,1)=x;
                epi(n,2)=y;
                n = n+1;
 end
 
 endoT = endo;
 epiT = epi;
 close(h) 