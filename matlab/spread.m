function [diout]=spread(idata,qdata,code1)
% 
% spread code sequence, similar to spread fortran function
% input : 
% idata, scalar or vector to be spread
% gdata, is the dimension 1 or 2
% code1, is the number of copies to replicate
% output:
% diout, the spread idata
%
% implemented by Diomar Cesar Lobao PhD
% UFF, Volta Redonda, RJ, Brazil
% 
switch nargin 
    case{0,1} 
        error('.....lack of input argument'); 
    case 2 
        code1=qdata; 
        qdata=idata; 
end 
%
[hn,vn]=size(idata) ;
% 
if (hn>vn || vn>hn)
    disp('array...'); 
    if(hn>vn)
        vc=hn;
    else
        vc=vn;
    end
else
    disp('scalar...');
    vc=hn;
end 
%
iout=zeros(vc,code1);
%
% spread the information correctly
 if (qdata == 2) % dim = 2
    for ii=1:code1
        iout(:,ii)=idata(1:vc); 
    end
    diout = iout;
 else
     if(qdata == 1) % dim = 1
         for ii=1:code1
             iout(:,ii)=idata(1:vc); 
         end
         diout = iout';
     end
 end
return 
% end of file. 
