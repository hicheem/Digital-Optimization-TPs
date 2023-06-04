
var -1..1:x1;
var -2..2:x2;
var -3..3:x3;

var int:y;
constraint y=(x1+x2+x3);

constraint x1+x2<2;
solve 
%:: int_search([x1,x2,x3], largest, indomain_max, complete)
:: int_search([x1,x2,x3], input_order, indomain_random, complete)
:: relax_and_reconstruct(array1d([x1,x2,x3]), 90)
:: restart_luby(10)
maximize y;
%varsel :  smallest; largest; dom_w_deg;input_order;  first_fail;

%valsel : indomain_min; indomain_max;  indomain_median;  indomain_random;