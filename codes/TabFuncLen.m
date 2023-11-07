function Table_Lenght =TabFuncLen( DX,G,S,ty )
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
   ii=0;
   X=0;
  if(ty~=2&&ty~=3)
    F=0.0;
    crtGain=G-0.005*G;
    while (F<crtGain)
      F=FUNC(X,G,S,ty);
      X=X+DX;
      ii=ii+1;
    end
   elseif (ty==3)
        ii=11;
   else
      crtGain=0.005*G;% Gaussian
       F=G;
  while (F>crtGain)
     ii=ii+1;
     F=FUNC(X,G,S,ty);
     X=X+DX;
  end
  end
 Table_Lenght=ii*2+1 ; % always odd
     
end

