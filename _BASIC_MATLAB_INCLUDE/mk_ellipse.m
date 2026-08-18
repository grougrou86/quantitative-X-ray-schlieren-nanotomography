
% mk_ellipse(xr,yr,x,y) => 2D ellipse (inner=0, outer=1);
% mk_ellipse(xr,yr,zr,x,y,z) => 3D ellipse (inner=0, outer=1);
% Revised by KR Lee 2012/6/27 , caring XR || YR < 0
% Revised by KR Lee 2012/6/29 , caring exact center of ellipse.

function H = mk_ellipse(varargin)
% mk_ellipse(xr,yr,x,y) => 2D ellipse (inner=0, outer=1);
% mk_ellipse(xr,yr,zr,x,y,z) => 3D ellipse (inner=0, outer=1);
Da =  cell2mat(varargin);

switch nargin
    case 4;
        XR=Da(1);
        YR=Da(2);
        X=Da(3);
        Y=Da(4);
        
        
        if XR > 0 &&  YR > 0
            [XX, YY]=meshgrid(1:X,1:Y);
            H = ((XX-(floor(X/2)+1))./XR).^2+((YY-(floor(Y/2)+1))./YR).^2>1.0;

        else 
            H = ones(Y,X);
           % H(round((X+1)/2),round((Y+1)/2)) =0;
        end

        return;
        
    case 6;
        XR=Da(1);
        YR=Da(2);
        ZR=Da(3);
        X=Da(4);
        Y=Da(5);
        Z=Da(6);
        
            if XR > 0 &&  YR > 0
            [XX, YY, ZZ]=meshgrid(1:X,1:Y,1:Z);
            H = ((XX-(floor(X/2)+1))./XR).^2+((YY-(floor(Y/2)+1))./YR).^2+((ZZ-(floor(Z/2)+1))./ZR).^2>1.0;

            else 
            H = ones(X,Y,Z);
            % H(round((X+1)/2),round((Y+1)/2)) =0;
            end
end