clear all

% read data
load('../data/polymodgb.csv')

% age cutoffs
n1 = 4;
n2 = 13;

b = zeros(3);

b(1,1) = sum(polymodgb(1:n1,1:n1),"all");
b(1,2) = sum(polymodgb(1:n1,n1+1:n2),"all");
b(1,3) = sum(polymodgb(1:n1,n2+1:end),"all");
b(2,1) = sum(polymodgb(n1+1:n2,1:n1),"all");
b(2,2) = sum(polymodgb(n1+1:n2,n1+1:n2),"all");
b(2,3) = sum(polymodgb(n1+1:n2,n2+1:end),"all");
b(3,1) = sum(polymodgb(n2+1:end,1:n1),"all");
b(3,2) = sum(polymodgb(n2+1:end,1:n1+1:n2),"all");
b(3,3) = sum(polymodgb(n2+1:end,n2+1:end),"all");

bsum = sum(b,"all");

contacts = b./bsum

save('./mats/contacts.mat',"contacts")