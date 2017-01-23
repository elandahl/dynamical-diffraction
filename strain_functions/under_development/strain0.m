function ss =strain0(s0,num_z,num_t)

% Wave headed away from surface is s1

for i = 1:num_z
    for j = 1:num_t
        if i-j < 0 && j-i < num_z
            s1(j,i) = -s0(1,j-i)/2;  % Reflection has a negative sign
        elseif i-j > 0 && i-j < num_z
            s1(j,i) = s0(1,i-j)/2;
        else
            s1(j,i) = 0;
        end
    end
end

% Wave headed towards surface is s2
% This wave undergoes reflection and flips sign at surface

for i = 1:num_z
    for j = 1:num_t
        if i+j > num_z
            s2(j,i) = 0; % No reflection on wave going into bulk
        elseif i+j < num_z
            s2(j,i) = s0(1,i+j)/2;
        else
            s2(j,i) = 0;
        end
    end
end

ss = s1+s2; % Add together positive and negative direction waves
end 