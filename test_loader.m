% Program to load all material properties as variables for a particular sample

load sample.dat;

sample_name = 'Si'; % Should be Si, Ge, GaAs, or InSb

ID = find(strcmp({sample.name}, sample_name)==1); % ID is code for sample
num_fields = length(fieldnames(sample(ID))); % Determine the number of fields
for i = 2:num_fields % Skip first field (it is always the name)
  vname=genvarname(getfield(getfield(sample(ID),fieldnames(sample(ID)){i}),'var'));
  vval = getfield(getfield(sample(ID),fieldnames(sample(ID)){i}),'val');
  eval([vname ' = vval;']);
end
