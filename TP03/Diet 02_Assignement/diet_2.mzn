% Diet problem. ____ STIGLER'S NUTRITION MODEL

%Parameters
int: n; %Number of nutrients
int: num_foods; %Number of foods
set of int: F; %set of foods
set of int: N; %set of nutrients
array[F, N] of int : a; %matrix: amount of nutrient j in food i, ∀i∈F, ∀j∈N
array[F] of int: ci; %cost per serving of food i, ∀i∈F
array[F] of int : Fmin;%minimum number of required servings of food i, ∀i∈F
array[F] of int :Fmax;%maximum allowable number of servings of food i, ∀i∈F
array[N] of int :Nmin;%minimum required level of nutrient j, ∀j∈N
array[N] of int :Nmax;%maximum allowable level of nutrient j, ∀j∈N
%Variables
%xi = number of servings of food i to purchase/consume, ∀i∈F
array[F] of var 0.0..11.0: xi;
var float: cost; constraint cost = sum(f in F)(ci[f]*xi[f]); 
%Objective Function:  the total cost of the food 

constraint cost>=0.0 /\ forall(nn in N) (
 (sum(f in F)(a[f,nn] * xi[f])) <= Nmax[nn] /\ (sum(f in F)(a[f,nn] * xi[f])) >= Nmin[nn] )/\                                                               
 forall(f in F) ( xi[f] <= Fmax[f] /\  xi[f] >= Fmin[f] );  
     
      
%Minimize
solve minimize cost;
%display
output[
  "number of servings of food f to be purchased daily:\n",
  "Riz : ", show(xi[1])," cost : $",show(xi[1]*ci[1]),"\n",
  "Lait : ", show(xi[2])," cost : $",show(xi[2]*ci[2]),"\n",
  "Pain : ", show(xi[3])," cost : $",show(xi[3]*ci[3]),"\n",
  "quantite des nutritions:\n",
  "vitamine A : ",show(a[1,1]*xi[1]+a[2,1]*xi[2]+a[3,1]*xi[3]),"\n",
  "Calories : ",show(a[1,2]*xi[1]+a[2,2]*xi[2]+a[3,2]*xi[3]),"\n",
  "Total cost : $ ",show(cost),"\n"
];

