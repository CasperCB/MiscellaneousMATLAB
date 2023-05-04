% Model from Masetti et al, 1983
function [mobility] = mobilitySi_Masetti(carrierConc, field, ntype)
    field = field ./ 100;

    if ntype == true
        mob1 = 68.5; % cm^2/ V*s
        mob2 = 56.1; % cm^2/ V*s
        mobMax = 1414; % cm^2/ V*s
        C_r = 9.2 * 10^16; % cm^3
        C_s = 3.41 * 10^20; % cm^3
        alpha = 0.711;
        beta = 1.98;
        e_c = 1.95 * 10^4;
        betaE = 1;
        
        mobility = mob1 + ((mobMax - mob1)./ (1 + (carrierConc ./ C_r).^alpha)) - ((mob2)./ (1 + (C_s ./ carrierConc).^beta));

    elseif ntype == false
        mob1 = 44.9; % cm^2/ V*s
        mob2 = 29.0; % cm^2/ V*s
        mobMax = 470.5; % cm^2/ V*s
        C_r = 2.23 * 10^17; % cm^3
        C_s = 6.10 * 10^20; % cm^3
        alpha = 0.719;
        beta = 2;
        p_c = 9.23 * 10^16; % cm^3
        e_c = 8 * 10^3;
        betaE = 2;
        
        mobility = mob1.*exp(-p_c ./ carrierConc) + ((mobMax)./ (1 + (carrierConc ./ C_r).^alpha)) - ((mob2)./ (1 + (C_s ./ carrierConc).^beta));
    
    end

    mobility0 = mobility;

    for i = 1:size(field,1)
        if field(i) ~= 0
            mobility(:,i) = (mobility0) ./ (1 + (field(i) ./ e_c).^betaE).^(1 + 1/betaE);
        end
    end

end

%% Old model
%
% N-type
%         mobMin = 65;
%         mobMax = 1330;
%         Nref = 6.3 * 10^16;
%         alpha = 0.76;
%         e_c = 1.95 * 10^4;
%         beta = 1;
% 
%         mobility = ((mobMax - mobMin) ./ (1 + (carrierConc ./ Nref).^alpha)) + mobMin;
%
% P-type
%         mobMin = 47.7;
%         mobMax = 495;
%         Nref = 8.5 * 10^16;
%         alpha = 0.72;
%         e_c = 8 * 10^3;
%         beta = 2;
% 
%         mobility = ((mobMax - mobMin) ./ (1 + (carrierConc ./ Nref).^alpha)) + mobMin;