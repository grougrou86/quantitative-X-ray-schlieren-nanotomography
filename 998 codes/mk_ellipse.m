%  H = mk_ellipse(XR,YR,XX,YY)
% Revised by KR Lee 2012/6/27 , caring XR || YR < 0
% Revised by KR Lee 2012/6/29 , caring exact center of ellipse.
% Revised by KR Lee 2018/1/10 , memory saving (2D only)

function H = mk_ellipse(varargin)
switch nargin
    case 4
        XR=varargin{1};
        YR=varargin{2};
        xx=varargin{3};
        yy=varargin{4};
        
        
        if (XR > 0) &&  (YR > 0)
            
            yvec = (1:1:yy)-(floor(yy/2)+1);
            xvec = (1:1:xx)-(floor(xx/2)+1);            
            
            
            [XX,YY]=meshgrid(xvec/XR,yvec/YR);
            [~,rho] = cart2pol(XX,YY);            
            H = rho>1.0;
%             
%             H = ((XX-(floor(xx/2)+1))./XR).^2+((YY-(floor(yy/2)+1))./YR).^2+((ZZ-(floor(zz/2)+1))./ZR).^2>1.0;
% 
%             
%             
%             coor =uint64([]);
%             cy = floor(YY/2)+1;
%             cx = floor(XX/2)+1;
%             for yd = 0:1:min(floor(YR),floor(YY/2))
%                 xd = floor(XR*sqrt(1-yd^2/YR^2));
%                 
%                 yvec = [yd*ones(1,2*xd+1),-yd*ones(1,2*xd+1)]+cy;
%                 xvec = [-xd:1:xd,-xd:1:xd]+cx;                
%                 validvec = logical((yvec > 0) .* (yvec <= YY) .* (xvec > 0) .* (xvec <= XX));                
%                 yvec = yvec(validvec);
%                 xvec = xvec(validvec);
%                 
%                 coor = [coor,sub2ind([YY,XX],yvec,xvec)];
%             end            
%             H = true(YY,XX);
%             H(coor)=0;                             
        else 
            H = true(yy,xx);
        end
        
    case 6
        XR=varargin{1};
        YR=varargin{2};
        ZR=varargin{3};
        xx=varargin{4};
        yy=varargin{5};
        zz=varargin{6};
        
            if XR > 0 &&  YR > 0
            [XX, YY, ZZ]=meshgrid(1:xx,1:yy,1:zz);
            H = ((XX-(floor(xx/2)+1))./XR).^2+((YY-(floor(yy/2)+1))./YR).^2+((ZZ-(floor(zz/2)+1))./ZR).^2>1.0;

            else 
            H = true(yy,xx,zz);
            % H(round((X+1)/2),round((Y+1)/2)) =0;
            end
end