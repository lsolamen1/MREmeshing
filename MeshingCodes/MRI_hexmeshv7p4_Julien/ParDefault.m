function [usedef]=ParDefault(nlabs,par,NbrInputs,default)
disp('MRE-Zone v7.3 Data Converter')
if(par)
    disp(['Using Parallelized version with ' int2str(nlabs) ' labs'])
else
    disp('Using non-Parallelized version - modify value of ''par'' to change')
end

if(NbrInputs==0)       % Default value is to prompt for inputs.
    default='no';
end
usedef=strcmp(default,'default'); % If no => false;else=>true
if(usedef)
    disp('Using Default values, no input prompts')
end
end
