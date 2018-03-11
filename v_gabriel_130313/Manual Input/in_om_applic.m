function [om_appl]= in_om_applic() 

value = 750;
value2 = value/10;

om_appl=[0 0/10^5 0/10^5
   168 value/10^5 value2/10^5
    534 value/10^5 value2/10^5
    890 value/10^5 value2/10^5];

%om_appl =[0 0 0];