function value=FUNC(X,G,S,ty)
%UNTITLED4 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
switch ty
    case -1
      value=2*G/(1+exp(-X*S))-G; % Signed Sigmoid
    case 0
       value=X*S;% Linear
    case 1
       value=G/(1+exp(-X*S));% Unsigned Sigmoidal
    case 2
       value=G*exp(-(X*X)/5.0); % Gaussian
    case 3
       value=2*G*(rand-0.5); % Random
    otherwise
       value=0;
end

